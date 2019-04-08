import time
import threading
import obspy
import httplib2
import logging
import numpy as np
import tensorflow as tf
import pymysql
import csv
import matplotlib.pyplot as plt

from config.config import Config
from eqdetector.eqdetector import EqDetector
from eqpickup.eqpickup import EqPickup
from location.location import Hyposat


class RequestData(object):
    '''
    从AWS中获取定时获取数据。
    '''

    def __init__(self):
        self.config = Config()
        self.stream_buffer = []
        self.data_buffer = []
        self.data_stats_buffer = []
        self.req_lock = threading.Lock()
        self.contain_windows = self.config.detect_window_size // self.config.window_lag_time - 1
        self.location_model = Hyposat()

        # 用于Cata_id中地震序号排序
        self.Cata_id_last_id = ''

    def process_data(self, starttime_UTC, endtime_UTC):
        '将数据处理成np.array, 并记录状态'
        for i, stream in enumerate(self.stream_buffer):
            # for stream in self.stream_buffer:
            if not (stream[0].stats.starttime == stream[1].stats.starttime
                    and stream[0].stats.starttime == stream[2].stats.starttime):
                starttime = max([stream[0].stats.starttime,
                                 stream[1].stats.starttime,
                                 stream[2].stats.starttime])
                stream = stream.slice(starttime=starttime)

            if not (stream[0].stats.endtime == stream[1].stats.endtime
                    and stream[0].stats.endtime == stream[2].stats.endtime):
                endtime = min([stream[0].stats.endtime,
                               stream[1].stats.endtime,
                               stream[2].stats.endtime])
                stream = stream.slice(endtime=endtime)
            self.stream_buffer[i] = stream

            data_stats = {'id': stream[0].id[:-4], 'starttime': stream[0].stats.starttime,
                          'endtime': stream[0].stats.endtime}
            self.data_stats_buffer.append(data_stats)

            data = np.array([stream[0].data, stream[1].data, stream[2].data])
            start_supplement_num = round((data_stats['starttime'] - starttime_UTC) * 100)
            end_supplement_num = round((endtime_UTC - data_stats['endtime']) * 100)

            if start_supplement_num != 0:
                start_supplement_data = np.zeros([3, start_supplement_num]) + np.mean(data[:, :200], axis=1,
                                                                                      keepdims=True)
                data = np.concatenate((start_supplement_data, data), axis=1)

            if end_supplement_num != 0:
                end_supplement_data = np.zeros([3, end_supplement_num]) + np.mean(data[:, -200:], axis=1, keepdims=True)
                data = np.concatenate((data, end_supplement_data), axis=1)
            self.data_buffer.append(data.T)

        self.data_buffer = np.array(self.data_buffer)

    def request_data(self, begin_time, end_time, is_plot=False, is_write2db=False):
        self.req_lock.acquire()

        print('Request data, time:{}'.format(begin_time))
        starttime_UTC, endtime_UTC = self.AWS_get_data(begin_time, end_time)
        print('Request data from AWS complete, time:{}'.format(
            time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))))
        if len(self.stream_buffer) != 0:
            self.process_data(starttime_UTC, endtime_UTC)
            self.analyse_data(starttime_UTC, endtime_UTC, is_plot, is_write2db)

        self.stream_buffer = []
        self.data_stats_buffer = []
        self.data_buffer = []

        self.req_lock.release()

    def AWS_get_data(self, begin_time, end_time):
        '从AWS中获取数据'
        h = httplib2.Http()
        host = self.config.request_data_config.jopens_host
        url = r'http://{}:8080/jopens-ws/app/aws/gwfw;'.format(host)
        starttime_UTC = obspy.UTCDateTime(begin_time) - 3600 * 8
        endtime_UTC = obspy.UTCDateTime(end_time) - 3600 * 8 - 0.01
        begin_time = begin_time.replace(' ', '%20')
        end_time = end_time.replace(' ', '%20')

        for station in self.config.request_data_config.monitor_station:
            net, sta, loc = station.split('.')

            tmp_url = url + 'net={};sta={};loc={};cha=.*;beg={};end={};'.format(net, sta, loc, begin_time, end_time)
            print(tmp_url)

            try:
                resp, content = h.request(tmp_url)

                with open('requestdata/aws_buffer_file', 'wb') as f:
                    f.write(content)

                eq_stream = obspy.read('requestdata/aws_buffer_file')
                eq_stream.merge(fill_value='interpolate')
                eq_stream.sort()

                eq_stream = eq_stream.slice(starttime=starttime_UTC, endtime=endtime_UTC)
                eq_stream.detrend('spline', order=2, dspline=50)
                # eq_stream.plot()
                self.stream_buffer.append(eq_stream)

            except Exception as e:
                logger = logging.getLogger(__name__)
                logger.setLevel(level=logging.INFO)
                handler = logging.FileHandler("requestdata/requestdatalog.txt")
                handler.setLevel(logging.INFO)
                formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
                handler.setFormatter(formatter)
                logger.addHandler(handler)

                logger.info("Reading data from AWS error, station:{}".format(station))
                logger.exception(e)

        return starttime_UTC, endtime_UTC

    def get_station_info(self):
        host = self.config.request_data_config.jopens_host
        db = pymysql.connect(host=host, user='root', passwd='rootme', db='jopens')
        cur = db.cursor()
        cur.execute("select station, longitude, latitude, elevation from Blockette050")
        stations = [station for station in cur]
        cur.close()
        db.close()
        return stations

    def analyse_data(self, starttime_UTC, endtime_UTC, is_plot=False, is_write2db=False):
        detector = EqDetector()

        tf_config = tf.ConfigProto()
        tf_config.gpu_options.allow_growth = True
        sess = tf.Session(config=tf_config)
        detector.start_sess(sess)
        event_list, current_endtime_UTC = detector.detect(sess, self.data_buffer, self.data_stats_buffer,
                                                          starttime_UTC,
                                                          endtime_UTC)
        detector.close_sess(sess)
        print(event_list)
        self.contain_windows = min((endtime_UTC - current_endtime_UTC) // self.config.window_lag_time,
                                   self.config.detect_window_size // self.config.window_lag_time - 1)

        event_stats_list = []
        pickup_model = EqPickup()
        sess = tf.Session(config=tf_config)
        pickup_model.start_sess(sess)
        if len(event_list) != 0:
            station2index = dict()
            for index in range(len(self.data_stats_buffer)):
                station = self.data_stats_buffer[index]['id']
                station2index[station] = index

            for event in event_list:
                event_stream = []
                for station in event:
                    stream = self.stream_buffer[station2index[station[0]]].slice(starttime=station[1],
                                                                                 endtime=station[2])
                    event_stream.append(stream)
                event_stats = pickup_model.pickup(sess, event_stream)
                event_stats_list.append(event_stats)

                if is_plot:
                    self.plot_event(event_stream, event_stats)

        pickup_model.close_sess(sess)

        filter_event_stats_list = []
        for event_stats in event_stats_list:
            act_num = 0
            filter_event_stats = []
            for station_stats in event_stats:
                if station_stats['p_confidence'] >= 0.5:
                    act_num += 1
                    if station_stats['s_confidence'] <= 0.7:
                        del station_stats['s_confidence']
                        del station_stats['S']

                    filter_event_stats.append(station_stats)

            if act_num >= 3:
                filter_event_stats_list.append(filter_event_stats)

        event_stats_list = filter_event_stats_list

        event_solution_stats_list = []
        if len(event_stats_list) > 0:
            stations = self.get_station_info()
            self.location_model.create_stations_file(stations)

            for event_stats in event_stats_list:
                event_starttime_UTC = event_stats[0]['starttime']
                event_id = event_starttime_UTC.strftime('%Y%m%d%H%M%S')
                output_file_name = "{}".format(event_id)
                event_solution_stats = self.location_model.loc(event_stats, output_file_name)
                if event_solution_stats is not None:
                    event_solution_stats_list.append(event_solution_stats)
        print(event_solution_stats_list)
        if len(event_solution_stats_list) != 0:
            # 提取出分析出来的地震事件状态（[北京时间，纬度，经度，深度，震级，rms]）
            write_lines = []
            for event_solution_stats in event_solution_stats_list:
                event_source_time = (event_solution_stats['source_time_UTC'] + 8 * 3600).strftime(
                    '%Y-%m-%d %H:%M:%S')
                line = [event_source_time, event_solution_stats['lat'], event_solution_stats['lon'],
                        event_solution_stats['z'], event_solution_stats['magnitude'],
                        event_solution_stats['rms']]
                write_lines.append(line)

            # 将结果写入jopens数据库
            if is_write2db:
                host = self.config.request_data_config.jopens_host
                db = pymysql.connect(host=host, user='root', passwd='rootme', db='jopens')
                cur = db.cursor()
                for line in write_lines:
                    try:
                        write_Net_code = 'SC'
                        write_Save_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))
                        write_Operator = 'USTC'
                        write_O_time = line[0]
                        write_Epi_lat = line[1]
                        write_Epi_lon = line[2]
                        write_Epi_depth = line[3]

                        # Cata_id 排序
                        write_id_time = write_Net_code + '.' + time.strftime('%Y%m%d%H%M', time.strptime(line[0], '%Y-%m-%d %H:%M:%S'))
                        if write_id_time == self.Cata_id_last_id[:-5]:
                            write_id_time = write_id_time + '.%04d' % (int(self.Cata_id_last_id[-4:]) + 1)
                        else:
                            write_id_time = write_id_time + '.0001'
                        self.Cata_id_last_id = write_id_time
                        write_id = write_id_time + '.' + 'A.001'



                        if line[4] == '':
                            write_M = -99999
                        else:
                            write_M = float(line[4].split()[0]) - 1     # 需求：人为算小一级
                        write_Rms = line[5]

                        # 超过2.5级不往数据库里写
                        if write_M <= 2.5:
                            print(
                                "INSERT INTO Catalog(Auto_flag, id, Net_code, Save_time, Operator, O_time, Epi_lat, Epi_lon, Epi_depth, M, Rms) VALUES ('A', '%s', '%s', '%s', '%s', '%s', %s, %s, %s, %s, %s)" % (
                                    write_id, write_Net_code, write_Save_time, write_Operator, write_O_time, write_Epi_lat,
                                    write_Epi_lon, write_Epi_depth, write_M, write_Rms))
                            cur.execute(
                                "INSERT INTO Catalog(Auto_flag, id, Net_code, Save_time, Operator, O_time, Epi_lat, Epi_lon, Epi_depth, M, Rms) VALUES ('A', '%s', '%s', '%s', '%s', '%s', %s, %s, %s, %s, %s)" % (
                                    write_id, write_Net_code, write_Save_time, write_Operator, write_O_time, write_Epi_lat,
                                    write_Epi_lon, write_Epi_depth, write_M, write_Rms))
                            db.commit()
                    except:
                        print('写入数据库错误！')
                        db.rollback()

                db.close()

            # 将结果写入csv文件
            with open('event_record/events{}_{}.csv'.format((starttime_UTC + 8 * 3600).strftime('%Y%m%d%H%M%S'),
                                                            (endtime_UTC + 8 * 3600).strftime('%Y%m%d%H%M%S')), 'w',
                      newline='') as events_file:
                csv_writer = csv.writer(events_file)
                for line in write_lines:
                    csv_writer.writerow(line)

    def plot_event(self, event_stream, event_stats):
        subplot_num = 10
        data = []
        for stream in event_stream:
            data.append(stream[2].data)

        station_num = len(data)
        i = 0
        figure_num = 1
        while i < station_num:
            plt.figure(figure_num)
            current_plot_num = min([subplot_num, station_num - i])
            for j in range(current_plot_num):
                plt.subplot(current_plot_num, 1, j + 1)
                plt.plot(data[i + j], 'k')
                if event_stats[i + j]['p_confidence'] >= 0.5:
                    event_p_start = int((event_stats[i + j]['P'] - event_stats[i + j]['starttime']) * 100)
                    plt.axvline(event_p_start, color='b')
                    if event_stats[i + j]['s_confidence'] >= 0.7:
                        event_s_start = int((event_stats[i + j]['S'] - event_stats[i + j]['starttime']) * 100)
                        plt.axvline(event_s_start, color='g')
            i = i + subplot_num
            figure_num = figure_num + 1

        plt.show()

    def scanning_specified_time(self, begin_time, end_time, is_plot=False, is_write2db=False):
        begin_timestamp = time.mktime(time.strptime(begin_time, '%Y-%m-%d %H:%M:%S'))
        end_timestamp = time.mktime(time.strptime(end_time, '%Y-%m-%d %H:%M:%S'))
        current_timestamp = begin_timestamp
        while current_timestamp < end_timestamp:
            current_begin_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(current_timestamp))
            current_end_time = min(current_timestamp + 3600 * 24, end_timestamp)
            current_end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(current_end_time))
            self.request_data(current_begin_time, current_end_time, is_plot, is_write2db)
            current_timestamp += 3600 * 24

    def main_thread(self):
        ''

        while True:
            run_time = int(time.time())
            begin_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(
                run_time - self.config.request_data_config.buffer_time - self.config.request_data_config.buffer_window_time - self.contain_windows * self.config.window_lag_time))
            end_time = time.strftime('%Y-%m-%d %H:%M:%S',
                                     time.localtime(run_time - self.config.request_data_config.buffer_time))

            get_data_thread = threading.Thread(target=self.request_data, kwargs={'begin_time': begin_time,
                                                                                 'end_time': end_time,
                                                                                 'is_write2db': True})
            get_data_thread.start()
            time.sleep(self.config.request_data_config.request_time)


if __name__ == '__main__':
    sess = tf.Session()
    model = EqDetector()
    model.start_sess(sess)
    model.close_sess(sess)
    tf.reset_default_graph()
    model = EqPickup()
    sess = tf.Session()
    model.start_sess(sess)
