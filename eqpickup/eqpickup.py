import numpy as np
import tensorflow as tf

from eqpickup.dl_tradition_model import DL_Tradition_Model

class EqPickup(object):
    def __init__(self):
        self.pickup_model = DL_Tradition_Model()

    def start_sess(self, sess):
        self.pickup_model.start_sess(sess)

    def close_sess(self, sess):
        self.pickup_model.close_sess(sess)


    def pickup(self, sess, event_stream):

        event_data = []
        data_stats = []
        for station in event_stream:
            station_data = np.array([station[0].data, station[1].data, station[2].data]).T
            event_data.append(station_data)
            station_stats = {'id': station[0].id[:-4], 'starttime': station[0].stats.starttime, 'endtime': station[0].stats.endtime}
            data_stats.append(station_stats)
        event_data = np.array(event_data)
        event_stats = self.pickup_model.pickup(sess, event_data, data_stats)
        return event_stats




