
class RequestDataConfig(object):
    def __init__(self):
        self.jopens_host = '10.51.141.18'
        self.request_time = 3600 * 1   # 请求数据时间间隔，每隔多久请求一次数据
        self.monitor_station_txt = 'config/monitor_station.txt'
        self.monitor_station = []
        with open(self.monitor_station_txt) as f:
            for station in f.readlines():
                station = station.strip('\n')
                self.monitor_station.append(station)
        # self.monitor_station = self.monitor_station[20:23]


        self.buffer_time = 3600 * 24  # 缓冲时间，处理数据与当前时间的差
        self.buffer_window_time = self.request_time  # 缓冲区时间长度，处理数据时间长度


class EqDetectorConfig(object):
    def __init__(self):
        self.max_window_size = 60

class DLTraditionModelConfig(object):
    def __init__(self):
        self.cblstm_step_size = 100
        self.cblstm_lstm_layer_num = 2
        self.cblstm_max_grad_norm = 5
        self.cblstm_keep_prob = 0.5
        self.cblstm_class_num = 3
        self.ps_lwinsize = 300      # s * 100
        self.ps_rwinsize = 100      # s * 100
        self.min_p_s_time = 25      # s * 100


class Config(object):
    def __init__(self):
        self.request_data_config = RequestDataConfig()
        self.eq_detector_config = EqDetectorConfig()
        self.dl_tradition_model_config = DLTraditionModelConfig()
        self.window_lag_time = 15   # 扫描窗口时间跨度
        self.detect_window_size = 45
        self.pickup_window_size = 60






if __name__ == '__main__':
    config = Config()




