import numpy as np

from obspy.signal.filter import bandpass
from config.config import Config
from eqdetector.cnn import CNN


def window_slide(length, win_len, win_step):
    i = length - win_len
    for i in range(0, length - win_len + 1, win_step):
        yield i
    if i + win_len != length:
        yield length - win_len


class EqDetector(object):
    def __init__(self):
        self.config = Config()
        self.detect_model = CNN()

    def start_sess(self, sess):
        self.detect_model.start_sess(sess)

    def close_sess(self, sess):
        self.detect_model.close_sess(sess)

    def data_preprocess(self, data, filter_model='bandpass', isabs=True):
        if filter_model == 'bandpass':
            data = bandpass(data, freqmin=4, freqmax=49, df=100)

        data = (data - np.mean(data)) / (np.max(np.absolute(data)) + 1)
        if isabs:
            data = np.absolute(data)
        return data

    def detect(self, sess, eq_data, eq_stats, starttime_UTC, endtime_UTC):
        first_start = 0
        act_flag = 0
        event_list = []
        station_act_flag = np.zeros([len(eq_stats)], dtype=np.int)  # 0:未检测到地震信号；1：当前检测的地震信号；-1：已检测到地震信号且已结束
        station_confidence = np.zeros([len(eq_stats)], dtype=np.float32)

        current_endtime = -1
        print(eq_data.shape)
        window_length = eq_data.shape[-2]
        window_size_sample_num = self.config.detect_window_size * 100
        window_lag_sample_num = self.config.window_lag_time * 100
        max_window_sample_num = self.config.eq_detector_config.max_window_size * 100

        for window_start in window_slide(window_length, window_size_sample_num, window_lag_sample_num):
            data_input = np.array(eq_data[:, window_start:window_start + window_size_sample_num, :])
            class_pred, confidence = self.detect_model.classify(sess=sess, input_=data_input)

            if window_start < current_endtime:
                continue

            if act_flag == 0 and sum(class_pred) < 3:
                continue

            if act_flag == 0 and sum(class_pred) >= 3:  # fisrt act
                first_start = window_start
                act_flag = 1
                for station in range(len(class_pred)):
                    if class_pred[station] == 1:
                        station_act_flag[station] = 1
                        if confidence[station] > station_confidence[station]:
                            station_confidence[station] = confidence[station]
            elif act_flag == 1:
                stop_flag = 1
                for station in range(len(class_pred)):
                    if station_act_flag[station] == 0 and class_pred[station] == 1 and window_start + window_size_sample_num - first_start <= max_window_sample_num:
                        station_act_flag[station] = 1
                        stop_flag = 0
                        if confidence[station] > station_confidence[station]:
                            station_confidence[station] = confidence[station]
                    elif station_act_flag[station] == 1:
                        if class_pred[station] == 0:
                            station_act_flag[station] = -1
                        elif class_pred[station] == 1:
                            stop_flag = 0
                            if confidence[station] > station_confidence[station]:
                                station_confidence[station] = confidence[station]
                if stop_flag == 1 or window_start + window_size_sample_num - first_start >= max_window_sample_num:
                    event = []
                    for station in range(len(class_pred)):
                        if station_act_flag[station] != 0:
                            station_record = [eq_stats[station]['id'],
                                              starttime_UTC + first_start / 100.0,
                                              min(starttime_UTC + (
                                                      window_start + self.config.pickup_window_size * 100 - 1) / 100.0,
                                                  endtime_UTC),
                                              station_confidence[station]]
                            event.append(station_record)
                    event_list.append(event)

                    station_act_flag = np.zeros([len(eq_stats)], dtype=np.int)
                    act_flag = 0
                    first_start = 0
                    station_confidence = np.zeros([len(eq_stats)], dtype=np.float32)
                    current_endtime = window_start + window_size_sample_num - 1

        current_endtime_UTC = starttime_UTC + current_endtime / 100.0
        return event_list, current_endtime_UTC
