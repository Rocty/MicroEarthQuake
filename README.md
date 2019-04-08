# MicroEarthQuake

#### 介绍
**三峡水库微震检测项目**
基于TensorFlow深度学习框架的微震检测系统

#### 安装教程
Ⅰ.无GPU版本：（推荐）
1.安装docker:
参考https://docs.docker.com/install/linux/docker-ce/ubuntu/#set-up-the-repository
和http://www.cnblogs.com/youran-he/p/9566059.html

```
$ sudo apt-get update
$ sudo apt-get install \
    apt-transport-https \
    ca-certificates \
    curl \
software-properties-common

$ sudo cp /etc/apt/sources.list /etc/apt/sources.list.bak
$ sudo sed -i 's/archive.ubuntu.com/mirrors.ustc.edu.cn' /etc/apt/sources.list
$ sudo apt update
$ curl -fsSL https://mirrors.ustc.edu.cn/docker-ce/linux/ubuntu/gpg | sudo apt-key add -
$ sudo add-apt-repository "deb [arch=amd64] https://mirrors.ustc.edu.cn/docker-ce/linux/ubuntu $(lsb_release -cs) stable"
$ sudo apt update
$ sudo apt install docker-ce
```

2.拉取镜像：

```
$ sudo docker login
```
账号：ustclinjf
密码：ustclinjf

若CPU较新（支持AVX2.0指令集）：
`$ sudo docker pull ustclinjf/microearthquake:v0.1-cpu`
若CPU较旧：
`$ sudo docker pull ustclinjf/microearthquake:v0.1-cpu-tf1.5`

3.修改台站数据信息：
具体修改代码文件MicroEarthQuake/config/monitor_station.txt文件，将其修改为需监控的台站名，格式参照该文件。（注意：目前只监控00台站，所以只填写00台站）

再修改代码文件MicroEarthQuake/config/config.py中
self.jonpens_host = ‘10.51.141.18’
后面IP地址改为jopens服务器地址。

4.修改镜像时间：
1)先运行容器：
$ sudo docker run -it --net=host -v /home/ahdy/Downloads/MicroEarthQuakce:/root/code --privileged=true ustclinjf/microearthquake:v0.1-cpu bash
其中/home/ahdy/Downloads/MicroEarthQuake为存放代码文件地址，ustclinjf/microearthquake:v0.1-cpu为第二步容器版本，需根据实际情况更改。
2)将容器放置后台运行：按Ctrl+P+Q（在Teamviewer上可能有按键冲突，需重复按几次才能生效）。
3)查看容器ID：
$ sudo docker ps –a
记下IMAGE为ustclinjf/microearthquake:{版本号}，STATUS距离当前最近的CONTAINER ID。
4)修改时间：
$ sudo docker cp /etc/localtime {containerid}:/etc/localtime

5.修改hyposat权限：
进入代码文件MicroEarthQuake/location/hyposat.6_0d/bin
然后执行：
$ sudo chmod 777 hyposat

6.重新进入docker容器：
$ sudo docker attach {记录的容器ID}

7.执行程序：
执行代码文件MicroEarthQuake/main.py，已挂在在容器环境中的/root/code：
$ python /root/code/main.py


Ⅱ.有GPU版本：
1.安装docker:同上
2.安装显卡驱动：https://blog.csdn.net/qq_37935670/article/details/80377196
3.安装nvidia-docker
4.拉取镜像：
$ sudo docker login
账号：ustclinjf
密码：ustclinjf
$ sudo docker pull ustclinjf/microearthquake:v0.1
5.修改镜像时间：
将容器版本修改为ustclinjf/microearthquake:v0.1
运行容器命令：
$ sudo docker run -it --runtime=nvidia --net=host -v /home/ahdy/Downloads/MicroEarthQuakce:/root/code --privileged=true ustclinjf/microearthquake:v0.1 bash

6修改hyposat权限：同上

7.重新进入docker容器：同上

8.执行程序：同上



#### 启动流程
以下每行以 $ 开头，即 $ 开始为一行，中间不断行。{XXX}需根据不同机器环境做修改（{}花括号不用打）。

1.查看容器ID
$ sudo docker ps -a
记下IMAGE为ustclinjf/microearthquake:{版本号}，STATUS距离当前最近的CONTAINER ID。

2.启动docker 环境：
$ sudo docker attach {上一步记录的ID}
（若容器已经停止运行，即status不是up，则需先执行：
$ sudo docker start {上一步记录的ID}
）

3.执行程序：
$ python /root/code/main.py
(若找不到该文件，改为: $ python /root/code/MicroEarthQuake/main.py)

注：版本号视机器而定：
1.有CUDA® 计算能力为 3.5 或更高的 NVIDIA® GPU 卡，为v0.1
2.CPU较新（支持AVX2.0指令集），为v0.1-cpu
3.CPU较旧，为v0.1-cpu-tf1.5