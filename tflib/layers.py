import tensorflow as tf


def variable_summaries(var, name=None):
    """Attach a lot of summaries to a Tensor (for TensorBoard visualization)."""
    with tf.name_scope(name):
        with tf.name_scope('summaries'):
            mean = tf.reduce_mean(var)
            tf.summary.scalar('mean', mean)
            with tf.name_scope('stddev'):
                stddev = tf.sqrt(tf.reduce_mean(tf.square(var - mean)))
            std_summ = tf.summary.scalar('stddev', stddev)
            max_summ = tf.summary.scalar('max', tf.reduce_max(var))
            min_summ = tf.summary.scalar('min', tf.reduce_min(var))
            his_summ = tf.summary.histogram('histogram', var)

    tf.add_to_collection('train_summary', std_summ)
    tf.add_to_collection('train_summary', max_summ)
    tf.add_to_collection('train_summary', min_summ)
    tf.add_to_collection('train_summary', his_summ)


def conv(input,
         filter,
         strides,
         padding,
         acti_func=tf.nn.relu,
         wd=None,
         bias=None,
         name=None):
    with tf.variable_scope(name) as scope:
        # kernel = tf.get_variable('weight',
        #                          shape=filter,
        #                          dtype=tf.float32,
        #                          initializer=tf.truncated_normal_initializer(stddev=0.1))
        kernel = tf.Variable(initial_value=tf.truncated_normal(filter, stddev=0.1), name='weight')
        variable_summaries(kernel, 'weight')

        if wd is not None:
            weight_decay = tf.multiply(tf.nn.l2_loss(kernel), wd, name='weight_loss')
            tf.add_to_collection('losses', weight_decay)

        if bias is not None:
            # bias = tf.get_variable('bias',
            #                        filter[-1],
            #                        dtype=tf.float32,
            #                        initializer=tf.constant_initializer(bias))
            bias = tf.Variable(tf.constant(bias, shape=[filter[-1]]), name='bias')
            variable_summaries(bias, 'bias')

        convolution = tf.nn.conv2d(input, kernel, strides=strides, padding=padding)
        act = acti_func(convolution + bias, name='activation')
        his_summ = tf.summary.histogram('activations', act)
        tf.add_to_collection('train_summary', his_summ)
        return act


def pool(input,
         ksize=[1, 1, 3, 1],
         strides=[1, 1, 3, 1],
         padding='SAME',
         pool_func=tf.nn.max_pool,
         name=None):
    with tf.variable_scope(name) as scope:
        return pool_func(input, ksize=ksize, strides=strides, padding=padding, name=name)


def Unfold(input, name=None):
    with tf.variable_scope(name) as scope:
        num_batch, height, width, num_channels = input.get_shape()
        return tf.reshape(input, [-1, (height * width * num_channels).value])


def fc(input,
       output_dim,
       input_dim=None,
       acti_func=tf.nn.relu,
       wd=None,
       name=None):
    with tf.variable_scope(name) as scope:
        # input_dim = tf.shape(input)[1]
        if input_dim is None:
            num_batch, input_dim = input.get_shape()
            input_dim = input_dim.value
        weights = tf.get_variable('weight',
                                  shape=[input_dim, output_dim],
                                  dtype=tf.float32,
                                  initializer=tf.truncated_normal_initializer(stddev=0.1))
        # weights = tf.Variable(tf.truncated_normal(shape=[input_dim.value, output_dim], stddev=0.1), name='weight')
        variable_summaries(weights, 'weight')
        if wd is not None:
            weight_decay = tf.multiply(tf.nn.l2_loss(weights), wd, name='weight_loss')
            tf.add_to_collection('losses', weight_decay)

        bias = tf.get_variable('bias',
                               output_dim,
                               dtype=tf.float32,
                               initializer=tf.constant_initializer(0.0))
        # bias = tf.Variable(tf.constant(0.0, shape=[output_dim]), name='bias')
        variable_summaries(bias, 'bias')
        output = tf.matmul(input, weights) + bias
        output = acti_func(output)
        his_summ = tf.summary.histogram('activations', output)
        tf.add_to_collection('train_summary', his_summ)
        return output


def fold(input, shape, name=None):
    with tf.variable_scope(name) as scope:
        num_batch, input_dim = input.get_shape()
        return tf.reshape(input, [-1, shape[1].value, shape[2].value, shape[3].value])


def deconv(input,
           filter,
           output_shape,
           strides=[1, 1, 5, 1],
           padding='SAME',
           acti_func=tf.nn.relu,
           wd=None,
           bias=None,
           name=None):
    with tf.variable_scope(name) as scope:
        # kernel = tf.get_variable('weight',
        #                          shape=filter,
        #                          dtype=tf.float32,
        #                          initializer=tf.truncated_normal_initializer(stddev=0.1))
        kernel = tf.Variable(tf.truncated_normal(filter, stddev=0.1), name='weight')
        variable_summaries(kernel, 'weight')
        if wd is not None:
            weight_decay = tf.multiply(tf.nn.l2_loss(kernel), wd, name='weight_loss')
            tf.add_to_collection('losses', weight_decay)

        if bias is not None:
            # bias = tf.get_variable('bias',
            #                        filter[-2],
            #                        dtype=tf.float32,
            #                        initializer=tf.constant_initializer(bias))
            bias = tf.Variable(tf.constant(bias, shape=[filter[-2]]), name='bias')
            variable_summaries(bias, 'bias')

        convolution = tf.nn.conv2d_transpose(input, kernel, output_shape=output_shape, strides=strides, padding=padding)
        act = acti_func(convolution + bias)
        his_summ = tf.summary.histogram('activations', act)
        tf.add_to_collection('train_summary', his_summ)
        return act


def conv1d(input,
           filter,
           strides,
           padding,
           acti_func=tf.nn.relu,
           wd=None,
           bias=None,
           is_training=True,
           name=None):
    if not is_training:
        wd = None

    with tf.variable_scope(name) as scope:
        kernel = tf.get_variable('weight',
                                 shape=filter,
                                 dtype=tf.float32,
                                 initializer=tf.truncated_normal_initializer(stddev=0.1))
        if is_training:
            variable_summaries(kernel, 'weight')

        if wd is not None:
            weight_decay = tf.multiply(tf.nn.l2_loss(kernel), wd, name='weight_loss')
            tf.add_to_collection('losses', weight_decay)

        if bias is not None:
            bias = tf.get_variable('bias',
                                   filter[-1],
                                   dtype=tf.float32,
                                   initializer=tf.constant_initializer(bias))
            if is_training:
                variable_summaries(bias, 'bias')

        convolution = tf.nn.conv1d(input, kernel, stride=strides, padding=padding)
        act = acti_func(convolution + bias, name='activation')
        if is_training:
            his_summ = tf.summary.histogram('activations', act)
            tf.add_to_collection('train_summary', his_summ)
        return act


def pool1d(input,
           ksize,
           strides,
           padding='SAME',
           pooling_type='MAX',
           name=None):
    return tf.nn.pool(input,
                      window_shape=ksize,
                      pooling_type=pooling_type,
                      padding=padding,
                      strides=strides,
                      name=name)


def Unfold1d(input, name=None):
    with tf.variable_scope(name) as scope:
        num_batch, width, num_channels = input.get_shape()
        return tf.reshape(input, [-1, (width * num_channels).value])


def birnn(input_, layers_num, num_units, keep_prob, seq_len, name=None):
    fw_cell = []
    bw_cell = []
    for i in range(layers_num):
        lstm_cell_fw = tf.nn.rnn_cell.LSTMCell(num_units=num_units,
                                               use_peepholes=True,
                                               forget_bias=1.0,
                                               state_is_tuple=True,
                                               reuse=tf.get_variable_scope().reuse)
        lstm_cell_fw = tf.nn.rnn_cell.DropoutWrapper(lstm_cell_fw, output_keep_prob=keep_prob)
        lstm_cell_bw = tf.nn.rnn_cell.LSTMCell(num_units=num_units,
                                               use_peepholes=True,
                                               forget_bias=1.0,
                                               state_is_tuple=True,
                                               reuse=tf.get_variable_scope().reuse)
        lstm_cell_bw = tf.nn.rnn_cell.DropoutWrapper(lstm_cell_bw, output_keep_prob=keep_prob)
        fw_cell.append(lstm_cell_fw)
        bw_cell.append(lstm_cell_bw)
    fw_cell = tf.nn.rnn_cell.MultiRNNCell(
        fw_cell, state_is_tuple=True)
    bw_cell = tf.nn.rnn_cell.MultiRNNCell(
        bw_cell, state_is_tuple=True)

    with tf.variable_scope(name):
        (outputs, _) = tf.nn.bidirectional_dynamic_rnn(cell_fw=fw_cell,
                                                       cell_bw=bw_cell,
                                                       inputs=input_,
                                                       sequence_length=seq_len,
                                                       dtype=tf.float32)
        output = tf.concat(outputs, 2)

    return output
