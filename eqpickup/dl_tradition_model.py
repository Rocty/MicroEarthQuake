import numpy as np
import copy
from spectrum import aryule
from obspy.signal.filter import bandpass
from phasepapy.phasepicker import aicdpicker
import math
import pywt
import matplotlib.pyplot as plt

from eqpickup.cblstm import CBLSTM
from config.config import Config


def data_process(data):
    def process(d):
        d = bandpass(d, freqmin=3, freqmax=15, df=100)

        d = (d - np.mean(d)) / (np.max(np.absolute(d)) + 0.00001)
        d = np.absolute(d)
        return d

    new_data = []
    for i in range(data.shape[0]):
        station_data = data[i].T
        station_data = np.array([process(d) for d in station_data]).T
        new_data.append(station_data)

    new_data = np.array(new_data)
    return new_data


def get_ar_aic(data, order=3, n_win=100, s_win=100):
    data_len = len(data)
    aic = np.zeros(data_len, dtype=np.float32)
    An = aryule(data[:n_win], order)[0]
    As = aryule(data[data_len - s_win:], order)[0]
    An = An * -1
    As = As * -1
    for i in range(data_len):
        if i <= order or i >= data_len - order:
            aic[i] = 1e7
            continue
        var_before = 0
        for j in range(order, i + 1):
            var_before += (data[j] - np.dot(An, data[j - order: j])) ** 2 / (i + 1 - order)
        var_after = 0
        for j in range(i + 1, data_len - order + 1):
            var_after += (data[j] - np.dot(As, data[j - order: j])) ** 2 / (data_len - order - i)
        aic[i] += (i - order) * np.log(var_before + 0.1) + (data_len - order - i) * np.log(var_after + 0.1)

    return aic


def get_energy(window):
    return np.sum(window ** 2)


def aic(data, stats, p_index, s_index, class_prob, ps_lwinsize, ps_rwinsize, min_p_s_time):
    try:
        data = data.T
        data_len = data.shape[1]
        data_raw = copy.deepcopy(data)

        bhz = bandpass(data[2], freqmin=3, freqmax=15, df=100)
        bhn = bandpass(data[1], freqmin=3, freqmax=15, df=100)
        bhe = bandpass(data[0], freqmin=3, freqmax=15, df=100)

        tmp_bhe = np.zeros([data_len], dtype=np.float32)
        tmp_bhn = np.zeros([data_len], dtype=np.float32)
        tmp_bhz = np.zeros([data_len], dtype=np.float32)

        tmp_bhe[0] = bhe[0]
        tmp_bhn[0] = bhn[0]
        tmp_bhz[0] = bhz[0]

        for i in range(data_len):
            if i == 0:
                continue
            tmp_bhz[i] = bhz[i] ** 2 + 3 * (bhz[i] - bhz[i - 1]) ** 2
            tmp_bhn[i] = bhn[i] ** 2 + 3 * (bhn[i] - bhn[i - 1]) ** 2
            tmp_bhe[i] = bhe[i] ** 2 + 3 * (bhe[i] - bhe[i - 1]) ** 2

        bhz = tmp_bhz
        bhn = tmp_bhn
        bhe = tmp_bhe

        p_range_start = 75
        p_range_end = data_len - 75 - 1
        p_cfd = -1
        if p_index != -1:
            p_range_start = max([p_index * 100 - ps_lwinsize, p_range_start])
            p_range_end = min([(p_index + 1) * 100 + ps_rwinsize, p_range_end])
            if p_index + 1 < len(class_prob):
                p_cfd = (class_prob[p_index, 1] + class_prob[p_index + 1, 1]) / 2.0
            else:
                p_cfd = class_prob[p_index, 1]

            p_win = []
            if sum(data_raw[2, p_range_start: p_range_end] == 0) <= 150:
                p_win.append(data_raw[2, p_range_start: p_range_end])
            elif sum(data_raw[1, p_range_start: p_range_end] == 0) <= 150:
                p_win.append(data_raw[1, p_range_start: p_range_end])
            elif sum(data_raw[0, p_range_start: p_range_end + 1] == 0) <= 150:
                p_win.append(data_raw[0, p_range_start: p_range_end + 1])
            if len(p_win) == 0:
                event_p_start = 0
                p_cfd = 0
            else:
                p_aic = get_ar_aic(p_win[0])
                event_p_start = p_range_start + np.argmin(p_aic)
        else:
            p_energy_ratio = np.zeros([data_len], dtype=np.float32) - 1.0
            i = p_range_start
            while i <= p_range_end:
                p_energy_ratio[i] = get_energy(bhz[i: i + 75 + 1]) / (get_energy(
                    bhz[i - 75: i + 1]) + 1e-12)
                i += 1
            event_p_start = np.argmax(p_energy_ratio)

        if p_cfd == -1:
            first = math.floor(event_p_start / 100)
            if first + 1 < len(class_prob):
                p_cfd = (class_prob[first, 1] + class_prob[first + 1, 1]) / 2.0
            else:
                p_cfd = class_prob[first, 1]

        s_range_start = event_p_start + min_p_s_time
        s_range_end = data_len - 75 - 1
        s_cfd = -1
        if s_index != -1:
            s_range_start = max([s_index * 100 - ps_lwinsize, s_range_start])
            s_range_end = min([(s_index + 1) * 100 + ps_rwinsize, s_range_end])
            if s_range_end <= s_range_start + 250:
                s_range_start = s_range_end - 251
                if s_range_start < 0:
                    s_range_start = 0
                    s_range_end = 251
            if s_index + 1 < len(class_prob):
                s_cfd = (class_prob[s_index, 2] + class_prob[s_index + 1, 2]) / 2.0
            else:
                s_cfd = class_prob[s_index, 2]

            s_bhe = data_raw[0, s_range_start: s_range_end]
            s_bhn = data_raw[1, s_range_start: s_range_end]
            s_bhz = data_raw[2, s_range_start: s_range_end]

            s_bhe = pywt.wavedec(s_bhe, 'db2', level=3)
            s_bhe = [s_bhe[0], None, None, None]
            s_bhe = pywt.waverec(s_bhe, 'db2')
            s_bhn = pywt.wavedec(s_bhn, 'db2', level=3)
            s_bhn = [s_bhn[0], None, None, None]
            s_bhn = pywt.waverec(s_bhn, 'db2')
            s_bhz = pywt.wavedec(s_bhz, 'db2', level=3)
            s_bhz = [s_bhz[0], None, None, None]
            s_bhz = pywt.waverec(s_bhz, 'db2')

            tmp_bhe = np.zeros([len(s_bhe)], dtype=np.float32)
            tmp_bhn = np.zeros([len(s_bhe)], dtype=np.float32)
            tmp_bhz = np.zeros([len(s_bhe)], dtype=np.float32)

            tmp_bhe[0] = s_bhe[0]
            tmp_bhn[0] = s_bhn[0]
            tmp_bhz[0] = s_bhz[0]

            for i in range(len(s_bhe)):
                if i == 0:
                    continue
                tmp_bhn[i] = s_bhn[i] ** 2 + 3 * (s_bhn[i] - s_bhn[i - 1]) ** 2
                tmp_bhe[i] = s_bhe[i] ** 2 + 3 * (s_bhe[i] - s_bhe[i - 1]) ** 2
                tmp_bhz[i] = s_bhz[i] ** 2 + 3 * (s_bhz[i] - s_bhz[i - 1]) ** 2

            s_bhn = tmp_bhn
            s_bhe = tmp_bhe
            s_bhz = tmp_bhz

            s_n_pick = 10000000
            s_e_pick = 10000000
            s_z_pick = 10000000
            if sum(s_bhn == 0) <= 150:
                s_n_aic = get_ar_aic(s_bhn, n_win=100, s_win=100)
                s_n_pick = np.argmin(s_n_aic)
            if sum(s_bhe == 0) <= 150:
                s_e_aic = get_ar_aic(s_bhe, n_win=100, s_win=100)
                s_e_pick = np.argmin(s_e_aic)
            if sum(s_bhz == 0) <= 150:
                s_z_aic = get_ar_aic(s_bhz, n_win=100, s_win=100)
                s_z_pick = np.argmin(s_z_aic)
            s_pick = np.min([s_e_pick, s_n_pick, s_z_pick])
            if s_pick == 10000000:
                event_s_start = 0
                s_cfd = 0
            else:
                event_s_start = s_range_start + s_pick
        else:
            s_energy_ratio = np.zeros([data_len], dtype=np.float32) - 1.0
            i = s_range_start
            while i <= s_range_end:
                bhe_energy_ratio = get_energy(bhe[i: i + 75] + 1) / (get_energy(
                    bhe[i - 75: i + 1]) + 1e-12)
                bhn_energy_ratio = get_energy(bhn[i: i + 75] + 1) / (get_energy(
                    bhn[i - 75: i + 1]) + 1e-12)
                s_energy_ratio[i] = (bhe_energy_ratio + bhn_energy_ratio) / 2.0
                i += 1
            event_s_start = np.argmax(s_energy_ratio)

        if s_cfd == -1:
            first = math.floor(event_s_start / 100)
            if first + 1 < len(class_prob):
                s_cfd = (class_prob[first, 2] + class_prob[first + 1, 2]) / 2.0
            else:
                s_cfd = class_prob[first, 2]

        if event_p_start + 300 + 1 >= len(bhz):
            event_p_start = len(bhz) - 300 - 2
        if event_p_start - 500 < 0:
            event_p_start = 500

        if event_s_start + 100 + 1 >= len(bhe):
            event_s_start = len(bhe) - 100 - 2
        if event_s_start - 100 < 0:
            event_s_start = 100

        ampl_n = max(abs(data_raw[1]))
        ampl_e = max(abs(data_raw[0]))
        ampl = max(ampl_n, ampl_e) / 1000.0

        p_snr = np.sqrt(get_energy(data_raw[2][event_p_start: event_p_start + 300 + 1])) / len(
            data_raw[2][event_p_start: event_p_start + 300 + 1]) / (np.sqrt(get_energy(
            data_raw[2][event_p_start - 300: event_p_start + 1]) + 1e-12) / len(
            data_raw[2][event_p_start - 300: event_p_start + 1]))

        s_snr = np.sqrt(get_energy(data_raw[2][event_p_start: event_p_start + 200 + 1])) / len(
            data_raw[2][event_p_start: event_p_start + 200 + 1]) / (np.sqrt(get_energy(
            data_raw[2][event_p_start - 200: event_p_start + 1]) + 1e-12) / len(
            data_raw[2][event_p_start - 200: event_p_start + 1]))

        # plt.subplot(4, 1, 1)
        # plt.plot(data_raw[0], 'k')
        # plt.axvline(event_p_start, color='b')
        # plt.axvline(event_s_start, color='g')
        # plt.legend(['BHE', 'P-wave', 'S-wave'], loc=2)
        # # plt.axvline(p_range_start, color='r', ls='--')
        # # plt.axvline(p_range_end, color='r', ls='--')
        # # plt.axvline(s_range_start, color='r', ls='--')
        # # plt.axvline(s_range_end, color='r', ls='--')
        # # plt.subplot(7, 1, 2)
        # # plt.plot(bhe, 'k')
        # # plt.axvline(event_p_start, color='b')
        # # plt.axvline(event_s_start, color='g')
        # # plt.axvline(p_range_start, color='r', ls='--')
        # # plt.axvline(p_range_end, color='r', ls='--')
        # # plt.axvline(s_range_start, color='r', ls='--')
        # # plt.axvline(s_range_end, color='r', ls='--')
        # plt.subplot(4, 1, 2)
        # plt.plot(data_raw[1], 'k', label='BHN')
        # plt.axvline(event_p_start, color='b')
        # plt.axvline(event_s_start, color='g')
        # plt.legend(['BHN'], loc=2)
        # # plt.axvline(p_range_start, color='r', ls='--')
        # # plt.axvline(p_range_end, color='r', ls='--')
        # # plt.axvline(s_range_start, color='r', ls='--')
        # # plt.axvline(s_range_end, color='r', ls='--')
        # # plt.subplot(7, 1, 4)
        # # plt.plot(bhn, 'k')
        # # plt.axvline(event_p_start, color='b')
        # # plt.axvline(event_s_start, color='g')
        # # plt.axvline(p_range_start, color='r', ls='--')
        # # plt.axvline(p_range_end, color='r', ls='--')
        # # plt.axvline(s_range_start, color='r', ls='--')
        # # plt.axvline(s_range_end, color='r', ls='--')
        # plt.subplot(4, 1, 3)
        # plt.plot(data_raw[2], 'k', label='BHZ')
        # plt.axvline(event_p_start, color='b')
        # plt.axvline(event_s_start, color='g')
        # plt.legend(['BHZ'], loc=2)
        # # plt.axvline(p_range_start, color='r', ls='--')
        # # plt.axvline(p_range_end, color='r', ls='--')
        # # plt.axvline(s_range_start, color='r', ls='--')
        # # plt.axvline(s_range_end, color='r', ls='--')
        # # plt.subplot(7, 1, 6)
        # # plt.plot(bhz, 'k')
        # # plt.axvline(event_p_start, color='b')
        # # plt.axvline(event_s_start, color='g')
        # # plt.axvline(p_range_start, color='r', ls='--')
        # # plt.axvline(p_range_end, color='r', ls='--')
        # # plt.axvline(s_range_start, color='r', ls='--')
        # # plt.axvline(s_range_end, color='r', ls='--')
        # plt.subplot(4, 1, 4)
        # plt.plot(class_prob[: int(data_len // 100) + 1, 0], 'b')
        # plt.plot(class_prob[: int(data_len // 100) + 1, 1], 'r')
        # plt.plot(class_prob[: int(data_len // 100) + 1, 2], 'k')
        # plt.axvline(p_index, color='r', ls='--')
        # plt.axvline(s_index, color='k', ls='--')
        # plt.legend(['noise-prob', 'P-prob', 'S-prob'], loc=2)
        # plt.show()

        p_time = stats['starttime'] + event_p_start / 100.0
        s_time = stats['starttime'] + event_s_start / 100.0
        stats['P'] = p_time
        stats['S'] = s_time
        stats['p_confidence'] = p_cfd
        stats['s_confidence'] = s_cfd
        stats['amplitude'] = ampl
        stats['p_snr'] = p_snr
        stats['s_snr'] = s_snr
    except Exception as e:
        print('error!!!')
        print(('error type: {}'.format(e)))
        return []

    return stats


class DL_Tradition_Model(object):
    def __init__(self):
        self.config = Config()
        self.DL_model = CBLSTM()

    def start_sess(self, sess):
        self.DL_model.start_sess(sess)

    def close_sess(self, sess):
        self.DL_model.close_sess(sess)

    def rough_pickup(self, sess, data):
        event_p_index, event_s_index, event_class_prob = self.DL_model.pickup_p_s(sess, data)
        return event_p_index, event_s_index, event_class_prob

    def precise_pickup(self, data, data_stats, event_p_index, event_s_index, event_class_prob):
        ps_lwinsize = self.config.dl_tradition_model_config.ps_lwinsize
        ps_rwinsize = self.config.dl_tradition_model_config.ps_rwinsize
        min_p_s_time = self.config.dl_tradition_model_config.min_p_s_time
        event_stats = []
        for i in range(len(data_stats)):
            event_stats.append(
                aic(data[i], data_stats[i], event_p_index[i], event_s_index[i], event_class_prob[i], ps_lwinsize,
                    ps_rwinsize, min_p_s_time))

        return event_stats

    def pickup(self, sess, data, data_stats):
        data_processed = data_process(data)
        event_p_index, event_s_index, event_class_prob = self.rough_pickup(sess, data_processed)

        event_stats = self.precise_pickup(data, data_stats, event_p_index, event_s_index, event_class_prob)
        return event_stats
