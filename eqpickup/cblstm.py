import tensorflow as tf
import numpy as np
import math
import sys
import os
from tensorflow.python import debug as tf_debug

from config.config import Config
from tflib.models import Model
from tflib import layers


def get_acc(pred_seq, targets, seq_len):
    acc_list = []
    for i, length in enumerate(seq_len):
        tmp_sum = sum(pred_seq[i, :length] == targets[i, :length])
        acc_list.append(tmp_sum * 1.0 / length)

    return sum(acc_list) * 1.0 / len(seq_len)


def get_arrival_index(sequence, max_len, wave_type):
    label = None
    if wave_type == 'P':
        label = 1
    elif wave_type == 'S':
        label = 2

    if label in sequence[: max_len]:
        return np.where(sequence[: max_len] == label)[0][0]
    else:
        return -1


def get_p_s_error(pred_seq, targets, seq_len):
    pred_p_index = []
    targets_p_index = []
    pred_s_index = []
    targets_s_index = []
    p_error = []
    s_error = []
    for i, length in enumerate(seq_len):
        tmp_pred_p_index = get_arrival_index(pred_seq[i], length, 'P')
        pred_p_index.append(tmp_pred_p_index)
        tmp_targets_p_index = get_arrival_index(targets[i], length, 'P')
        targets_p_index.append(tmp_targets_p_index)
        tmp_pred_s_index = get_arrival_index(pred_seq[i], length, 'S')
        pred_s_index.append(tmp_pred_s_index)
        tmp_targets_s_index = get_arrival_index(targets[i], length, 'S')
        targets_s_index.append(tmp_targets_s_index)

        if tmp_targets_p_index != -1:
            p_error.append(abs(tmp_targets_p_index - tmp_pred_p_index))
        if tmp_targets_s_index != -1:
            s_error.append(abs(tmp_targets_s_index - tmp_pred_s_index))

    p_err_ave = np.mean(p_error)
    p_err_max = max(p_error)
    if len(s_error) == 0:
        s_err_ave = 0
        s_err_max = 0
    else:
        s_err_ave = np.mean(s_error)
        s_err_max = max(s_error)
    return p_err_ave, p_err_max, s_err_ave, s_err_max


class CBLSTM(object):
    def __init__(self):
        self.config = Config()
        self.step_size = self.config.dl_tradition_model_config.cblstm_step_size
        self.reader = None
        self.layer = self.setup_layer()
        self.loss = self.setup_loss()
        self.train_op = self.setup_train_op()
        self.train_metrics = self.setup_metrics(True)
        self.test_metrics = self.setup_metrics(False)
        self.train_merged = tf.summary.merge(tf.get_collection('train_summary'))
        self.train_metrics_merged = tf.summary.merge(tf.get_collection('train_metrics_summary'))
        self.test_metrics_merged = tf.summary.merge(tf.get_collection('test_metrics_summary'))
        self.train_writer = tf.summary.FileWriter('eqpickup/summary/cblstm')

    def start_sess(self, sess):
        saver, global_step = Model.continue_previous_session(sess,
                                                             model_file='cblstm',
                                                             ckpt_file='eqpickup/saver/cblstm/checkpoint')

    def close_sess(self, sess):
        sess.close()
        tf.reset_default_graph()

    def data_padding_preprocess(self, data, data_type):
        step_size = self.step_size
        sequence_len = list(map(len, data))
        sequence_len = [math.ceil(i / float(step_size)) for i in sequence_len]
        max_len = max(sequence_len)
        if max_len % step_size != 0:
            max_len = math.ceil(max_len)

        if data_type == 'input':
            result = np.zeros([len(data), max_len * step_size, 3], dtype=np.float32)
            for i, example in enumerate(data):
                for j, row in enumerate(example):
                    for k, val in enumerate(row):
                        result[i][j][k] = val

        elif data_type == 'targets':
            result = np.zeros([len(data), max_len], dtype=np.int32)
            for i, example in enumerate(data):
                for step, val in enumerate(example[::step_size]):
                    result[i][step] = np.max(val)

        return result, sequence_len

    def setup_layer(self):
        # LSTM_units_num = self.config.cblstm_l
        lstm_layers_num = self.config.dl_tradition_model_config.cblstm_lstm_layer_num
        class_num = self.config.dl_tradition_model_config.cblstm_class_num
        # batch_size = self.config.cblstm_batch_size

        layer = dict()
        input_ = tf.placeholder(tf.float32, shape=[None, None, 3], name='input')
        sequence_length = tf.placeholder(tf.int32, shape=[None], name='seq_len')
        targets = tf.placeholder(tf.int32, shape=[None, None], name='targets')
        keep_prob = tf.placeholder(tf.float32, name='keep_prob')

        batch_size = tf.shape(input_)[0]

        layer['input'] = input_
        layer['seq_len'] = sequence_length
        layer['targets'] = targets
        layer['keep_prob'] = keep_prob

        layer['conv1'] = layers.conv1d(layer['input'],
                                       filter=[4, 3, 8],
                                       strides=1,
                                       padding='SAME',
                                       wd=5e-5,
                                       bias=0.0,
                                       name='conv1')
        layer['pooling1'] = layers.pool1d(layer['conv1'],
                                          ksize=[2],
                                          strides=[2],
                                          padding='SAME',
                                          name='pooling1')
        layer['conv2'] = layers.conv1d(layer['pooling1'],
                                       filter=[4, 8, 16],
                                       strides=1,
                                       padding='SAME',
                                       wd=5e-5,
                                       bias=0.0,
                                       name='conv2')
        layer['pooling2'] = layers.pool1d(layer['conv2'],
                                          ksize=[2],
                                          strides=[2],
                                          padding='SAME',
                                          name='pooling2')
        layer['unfold'] = tf.reshape(layer['pooling2'], [batch_size, -1, 400])
        layer['unfold'] = tf.reshape(layer['unfold'], [-1, 400])

        layer['unfold'] = tf.nn.dropout(layer['unfold'], keep_prob)

        layer['dim_red'] = layers.fc(layer['unfold'], output_dim=100, wd=5e-5, name='dim_red')

        layer['dim_red'] = tf.reshape(layer['dim_red'], [batch_size, -1, 100])

        lstm_cell_fw1 = tf.nn.rnn_cell.LSTMCell(num_units=100,
                                                forget_bias=1.0,
                                                state_is_tuple=True,
                                                reuse=tf.get_variable_scope().reuse)
        lstm_cell_fw2 = tf.nn.rnn_cell.LSTMCell(num_units=100,
                                                num_proj=50,
                                                forget_bias=1.0,
                                                state_is_tuple=True,
                                                reuse=tf.get_variable_scope().reuse)
        lstm_cell_bw1 = tf.nn.rnn_cell.LSTMCell(num_units=100,
                                                forget_bias=1.0,
                                                state_is_tuple=True,
                                                reuse=tf.get_variable_scope().reuse)
        lstm_cell_bw2 = tf.nn.rnn_cell.LSTMCell(num_units=100,
                                                num_proj=50,
                                                forget_bias=1.0,
                                                state_is_tuple=True,
                                                reuse=tf.get_variable_scope().reuse)
        lstm_cell_fw1 = tf.nn.rnn_cell.DropoutWrapper(lstm_cell_fw1, output_keep_prob=keep_prob)
        lstm_cell_fw2 = tf.nn.rnn_cell.DropoutWrapper(lstm_cell_fw2, output_keep_prob=keep_prob)
        lstm_cell_bw1 = tf.nn.rnn_cell.DropoutWrapper(lstm_cell_bw1, output_keep_prob=keep_prob)
        lstm_cell_bw2 = tf.nn.rnn_cell.DropoutWrapper(lstm_cell_bw2, output_keep_prob=keep_prob)

        cell_fw = tf.nn.rnn_cell.MultiRNNCell(
            [lstm_cell_fw1, lstm_cell_fw2], state_is_tuple=True)
        cell_bw = tf.nn.rnn_cell.MultiRNNCell(
            [lstm_cell_bw1, lstm_cell_bw2], state_is_tuple=True)

        layer['dim_red'] = tf.nn.dropout(layer['dim_red'], keep_prob)

        with tf.variable_scope('bi_rnn'):
            (outputs, _) = tf.nn.bidirectional_dynamic_rnn(cell_fw=cell_fw,
                                                           cell_bw=cell_bw,
                                                           inputs=layer['dim_red'],
                                                           sequence_length=sequence_length,
                                                           dtype=tf.float32)
            output = tf.concat(outputs, 2)
            layer['birnn'] = tf.reshape(output, [-1, 50 * 2])

        layer['birnn'] = tf.nn.dropout(layer['birnn'], keep_prob)
        layer['fc'] = layers.fc(layer['birnn'], output_dim=50, wd=5e-5, name='fc')

        with tf.variable_scope('softmax'):
            softmax_w = tf.get_variable(name='softmax_w',
                                        shape=[50, class_num],
                                        initializer=tf.truncated_normal_initializer(stddev=0.05),
                                        dtype=tf.float32)
            weight_decay = tf.multiply(tf.nn.l2_loss(softmax_w), 5e-5, name='weight_loss')
            tf.add_to_collection('losses', weight_decay)
            softmax_b = tf.get_variable(name='softmax_b',
                                        shape=[class_num],
                                        initializer=tf.constant_initializer(value=0),
                                        dtype=tf.float32)
            xw_plus_b = tf.nn.xw_plus_b(layer['fc'], softmax_w, softmax_b)
            logits = tf.reshape(xw_plus_b, [batch_size, -1, class_num])
            layer['logits'] = logits
            class_prob = tf.nn.softmax(logits)
            layer['class_prob'] = class_prob

        layer['pred_seq'] = tf.cast(tf.argmax(class_prob, axis=2), tf.int32)
        return layer

    def setup_loss(self):
        with tf.name_scope('loss'):
            # fake_loss = tf.nn.sparse_softmax_cross_entropy_with_logits(logits=self.training_layer['logits'],
            #                                                            labels=self.training_layer['targets'])
            # mask = tf.cast(tf.sign(self.training_layer['targets']), dtype=tf.float32)
            mask = tf.sequence_mask(self.layer['seq_len'], dtype=tf.float32)
            # loss_per_example_per_step = tf.multiply(fake_loss, mask)
            # loss_per_example_sum = tf.reduce_sum(loss_per_example_per_step, reduction_indices=[1])
            # loss_per_example_average = tf.div(x=loss_per_example_sum,
            #                                   y=tf.cast(self.training_layer['seq_len'], tf.float32))
            raw_loss = tf.contrib.seq2seq.sequence_loss(self.layer['logits'],
                                                        self.layer['targets'],
                                                        mask,
                                                        average_across_timesteps=True,
                                                        average_across_batch=True)
            # loss = tf.reduce_mean(loss_per_example_average, name='loss')
            # loss_per_example_per_step = tf.multiply(raw_loss, mask)
            # mask = tf.reduce_sum(mask, axis=0)
            # loss_per_example_per_step = tf.div(x=loss_per_example_per_step,
            #                                    y=mask)
            raw_loss_summ = tf.summary.scalar('raw_loss', raw_loss)
            tf.add_to_collection('losses', raw_loss)
            loss = tf.add_n(tf.get_collection('losses'), name='total_loss')
            loss_summ = tf.summary.scalar('total_loss', loss)
            tf.add_to_collection('train_summary', raw_loss_summ)
            tf.add_to_collection('train_summary', loss_summ)
        return loss

    def setup_train_op(self):
        tvars = tf.trainable_variables()
        grads, _ = tf.clip_by_global_norm(tf.gradients(self.loss, tvars),
                                          self.config.dl_tradition_model_config.cblstm_max_grad_norm)
        optimizer = tf.train.AdamOptimizer(5e-4)
        return optimizer.apply_gradients(zip(grads, tvars))
        # return optimizer.minimize(self.loss)

    def setup_metrics(self, is_training):
        if is_training:
            name_scope = 'training'
            summary_name = 'train_metrics_summary'
        else:
            name_scope = 'test'
            summary_name = 'test_metrics_summary'

        metrics = dict()
        with tf.name_scope(name_scope):
            with tf.variable_scope(name_scope):
                metrics['acc'] = tf.placeholder(dtype=tf.float32, name='acc')
                acc_summ = tf.summary.scalar('acc', metrics['acc'])
                metrics['p_error'] = tf.placeholder(dtype=tf.float32, name='p_error')
                p_err_summ = tf.summary.scalar('p_error', metrics['p_error'])
                metrics['p_error_max'] = tf.placeholder(dtype=tf.float32, name='p_error_max')
                p_err_max_summ = tf.summary.scalar('p_error_max', metrics['p_error_max'])
                metrics['s_error'] = tf.placeholder(dtype=tf.float32, name='s_error')
                s_err_summ = tf.summary.scalar('s_error', metrics['s_error'])
                metrics['s_error_max'] = tf.placeholder(dtype=tf.float32, name='s_error_max')
                s_err_max_summ = tf.summary.scalar('s_error_max', metrics['s_error_max'])
                tf.add_to_collection(summary_name, acc_summ)
                tf.add_to_collection(summary_name, p_err_summ)
                tf.add_to_collection(summary_name, p_err_max_summ)
                tf.add_to_collection(summary_name, s_err_summ)
                tf.add_to_collection(summary_name, s_err_max_summ)

        return metrics

    def train(self, passes, new_training=True):
        sess_config = tf.ConfigProto()
        sess_config.gpu_options.allow_growth = True
        with tf.Session(config=sess_config) as sess:
            if new_training:
                saver, global_step = Model.start_new_session(sess)
            else:
                saver, global_step = Model.continue_previous_session(sess,
                                                                     model_file='cblstm',
                                                                     ckpt_file='eqpickup/saver/cblstm/checkpoint')

            self.train_writer.add_graph(sess.graph, global_step=global_step)

            for step in range(1 + global_step, 1 + passes + global_step):
                with tf.variable_scope('Train'):
                    input_, targets = self.reader.get_cblstm_batch_data('train',
                                                                        self.reader.pre_data_generator,
                                                                        self.reader.pre_validation_data)
                    input_, seq_len = self.data_padding_preprocess(input_, 'input')
                    targets, _ = self.data_padding_preprocess(targets, 'targets')
                    _, train_summary, loss, pred_seq = sess.run(
                        [self.train_op, self.train_merged, self.loss, self.layer['pred_seq']],
                        feed_dict={self.layer['input']: input_,
                                   self.layer['targets']: targets,
                                   self.layer['seq_len']: seq_len,
                                   self.layer['keep_prob']: self.config.dl_tradition_model_config.cblstm_keep_prob})
                    self.train_writer.add_summary(train_summary, step)

                    train_p_err, train_p_err_max, train_s_err, train_s_err_max = get_p_s_error(pred_seq, targets,
                                                                                               seq_len)
                    train_acc = get_acc(pred_seq, targets, seq_len)

                    [train_metrics_summary] = sess.run(
                        [self.train_metrics_merged],
                        feed_dict={self.train_metrics['acc']: train_acc,
                                   self.train_metrics['p_error']: train_p_err,
                                   self.train_metrics['p_error_max']: train_p_err_max,
                                   self.train_metrics['s_error']: train_s_err,
                                   self.train_metrics['s_error_max']: train_s_err_max})
                    self.train_writer.add_summary(train_metrics_summary, step)
                    print("gobal_step {},"
                          " training_loss {},"
                          " accuracy {},"
                          " p_error {},"
                          " p_err_max {},"
                          " s_error {},"
                          " s_err_max {}.".format(step, loss, train_acc, train_p_err, train_p_err_max, train_s_err,
                                                  train_s_err_max))

                if step % 50 == 0:
                    with tf.variable_scope('Test', reuse=True):
                        test_input, test_targets = self.reader.get_cblstm_batch_data('test',
                                                                                     self.reader.pre_data_generator,
                                                                                     self.reader.pre_validation_data)
                        test_input, test_seq_len = self.data_padding_preprocess(test_input, 'input')
                        test_targets, _ = self.data_padding_preprocess(test_targets, 'targets')
                        [test_pred_seq] = sess.run([self.layer['pred_seq']],
                                                   feed_dict={self.layer['input']: test_input,
                                                              self.layer['seq_len']: test_seq_len,
                                                              self.layer['keep_prob']: 1.0})
                        test_p_err, test_p_err_max, test_s_err, test_s_err_max = get_p_s_error(test_pred_seq,
                                                                                               test_targets,
                                                                                               test_seq_len)
                        test_acc = get_acc(test_pred_seq,
                                           test_targets,
                                           test_seq_len)
                        [test_metrics_summary] = sess.run(
                            [self.test_metrics_merged],
                            feed_dict={self.test_metrics['acc']: test_acc,
                                       self.test_metrics['p_error']: test_p_err,
                                       self.test_metrics['p_error_max']: test_p_err_max,
                                       self.test_metrics['s_error']: test_s_err,
                                       self.test_metrics['s_error_max']: test_s_err_max})
                        self.train_writer.add_summary(test_metrics_summary, step)
                        print("test_acc {}, "
                              "test_p_err {},"
                              "test_p_err_max {},"
                              "test_s_err {},"
                              "test_s_err_max {}.".format(test_acc, test_p_err, test_p_err_max, test_s_err,
                                                          test_s_err_max))

                if step % 50 == 0:
                    saver.save(sess, 'eqpickup/saver/cblstm/cblstm', global_step=step)
                    print('checkpoint saved')

    def pickup_p_s(self, sess, input_, get_pred_seq=False):
        # with tf.Session() as sess:
        #     saver, global_step = Model.continue_previous_session(sess,
        #                                                          model_file='bi_rnn',
        #                                                          ckpt_file='saver/bi_rnn/checkpoint')
        with tf.variable_scope('model', reuse=True):
            input_, seq_len = self.data_padding_preprocess(input_, 'input')
            pred_seq, class_prob = sess.run([self.layer['pred_seq'], self.layer['class_prob']],
                                            feed_dict={self.layer['input']: input_,
                                                       self.layer['seq_len']: seq_len,
                                                       self.layer['keep_prob']: 1.0})
            p_index = [get_arrival_index(pred_seq[i], seq_len[i], 'P') for i in range(len(input_))]
            s_index = [get_arrival_index(pred_seq[i], seq_len[i], 'S') for i in range(len(input_))]
        if get_pred_seq:
            return p_index, s_index, class_prob, pred_seq
        else:
            return p_index, s_index, class_prob


if __name__ == '__main__':
    cblstm = CBLSTM()
    cblstm.train(13000, True)

    ### test acc
    # cblstm = CBLSTM()
    # sess_config = tf.ConfigProto()
    # sess_config.gpu_options.allow_growth = True
    # sess = tf.Session(config=sess_config)
    # saver, global_step = Model.continue_previous_session(sess,
    #                                                      model_file='cblstm',
    #                                                      ckpt_file='saver/cblstm/checkpoint')

    # acc_ave = 0
    # p_err_ave = 0
    # s_err_ave = 0
    # p_err_max = []
    # s_err_max = []
    # for i in range(20):
    #     test_input, test_targets = cblstm.reader.get_cblstm_batch_data('test')
    #     test_targets, test_seq_len = cblstm.data_padding_preprocess(test_targets, 'targets')
    #     events_p_index, events_s_index, _, pred_seq = cblstm.pickup_p_s(sess, test_input, get_pred_seq=True)
    #     tmp_acc = get_acc(pred_seq,
    #                       test_targets,
    #                       test_seq_len)
    #     acc_ave += tmp_acc
    #     tmp_p_error, tmp_p_err_max, tmp_s_error, tmp_s_err_max = get_p_s_error(pred_seq,
    #                                                                            test_targets,
    #                                                                            test_seq_len)
    #     p_err_ave += tmp_p_error
    #     s_err_ave += tmp_s_error
    #     p_err_max.append(tmp_p_err_max)
    #     s_err_max.append(tmp_s_err_max)
    #     print('acc:{}, p_err:{}, p_err_max:{}, s_err:{}, s_err_max:{}.'.format(tmp_acc, tmp_p_error, tmp_p_err_max,
    #                                                                            tmp_s_error, tmp_s_err_max))
    #
    # acc_ave /= 20
    # p_err_ave /= 20
    # s_err_ave /= 20
    # p_err_max = max(p_err_max)
    # s_err_max = max(s_err_max)
    # print('acc:{}, p_err:{}, p_err_max:{}, s_err:{}, s_err_max:{}.'.format(acc_ave, p_err_ave, p_err_max, s_err_ave,
    #                                                                        s_err_max))
