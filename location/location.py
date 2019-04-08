import os
import platform

from location.util import *

class Location(object):
    def __init__(self):
        pass

    def loc(self, event, output_filename):
        pass


class Hyposat(Location):
    def create_stations_file(self, stations):
        write_stations_dat(stations, r'location/hyposat.6_0d/data/stations.dat')

    def loc(self, event, output_filename):
        write_phases(event, 'HYPOSAT', r'location/hyposat.6_0d/hyposat-in')
        loc_output_filename = None
        if 'Windows' in platform.system():
            os.system("location\hyposat.6_0d\loc.bat")
            os.system("copy location\hyposat.6_0d\hyposat-out event_record\{}-out".format(output_filename))
            os.system("copy location\hyposat.6_0d\hyposat-isf event_record\{}-isf".format(output_filename))
            loc_output_filename = "event_record\{}-out".format(output_filename)
        elif 'Linux' in platform.system():
            os.system("bash location/hyposat.6_0d/run_ubuntu")
            os.system(r"cp location/hyposat.6_0d/hyposat-out event_record/{}-out".format(output_filename))
            os.system("cp location/hyposat.6_0d/hyposat-isf event_record/{}-isf".format(output_filename))
            loc_output_filename = "event_record/{}-out".format(output_filename)


        event_solution_stats = read_hyposat_out(loc_output_filename)
        return event_solution_stats




