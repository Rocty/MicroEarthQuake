import tensorflow as tf

from requestdata.requestdata import RequestData



if __name__ == '__main__':
    request_data = RequestData()
    request_data.main_thread()
    # request_data.scanning_specified_time('2019-01-14 05:00:00', '2019-01-14 08:00:00', is_plot=False)
    # request_data.scanning_specified_time('2019-01-14 15:41:00', '2019-01-14 15:44:00', is_plot=False)
    # request_data.scanning_specified_time('2019-01-14 16:50:00', '2019-01-14 16:54:00', is_plot=False)
    # request_data.scanning_specified_time('2019-01-23 00:00:00', '2019-01-24 00:00:00', is_plot=False)
