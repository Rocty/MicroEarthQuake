import math
import obspy


def write_phases(event, fformat, filename):
    if fformat == 'HYPOSAT':
        print("Writing phases to %s for HYPOSAT" % filename)
        with open(filename, 'w') as fid:
            net = event[0]['id'].split('.')[0]
            starttime = event[0]['starttime']
            endtime = event[0]['endtime']
            fid.write('net %s, event record starttime(UTC): %s, endtime(UTC): %s \n' % (net, starttime, endtime))
            for station in event:
                if 'P' in station.keys():
                    if station['p_confidence'] > 0.5:
                        station_name = station['id'].split('.')[1]
                        p_year = station['P'].year
                        p_month = station['P'].month
                        p_day = station['P'].day
                        p_hh = station['P'].hour
                        p_mm = station['P'].minute
                        p_ss = station['P'].second
                        p_ms = station['P'].microsecond
                        p_ss = p_ss + p_ms / 1000000.0
                        p_std = 0.1
                        amplitude = station['amplitude']
                        p_snr = station['p_snr']
                        fid.write(
                            '%-5s P1       %4.0f %02d %02d %02d %02d %05.2f  %5.3f -999.   0.00 -999.  0.00 T__DRM_ -999.  %12.2f %7.2f\n' % (
                                station_name, p_year, p_month, p_day, p_hh, p_mm, p_ss, p_std, amplitude, p_snr))

                if 'S' in station.keys():
                    if station['s_confidence'] > 0.5:
                        station_name = station['id'].split('.')[1]
                        s_year = station['S'].year
                        s_month = station['S'].month
                        s_day = station['S'].day
                        s_hh = station['S'].hour
                        s_mm = station['S'].minute
                        s_ss = station['S'].second
                        s_ms = station['S'].microsecond
                        s_ss = s_ss + s_ms / 1000000.0
                        s_std = 0.2
                        amplitude = station['amplitude']
                        s_snr = station['s_snr']
                        fid.write(
                            '%-5s S1       %4.0f %02d %02d %02d %02d %05.2f  %5.3f -999.   0.00 -999.  0.00 T__DRM_ -999.  %12.2f %7.2f\n' % (
                                station_name, s_year, s_month, s_day, s_hh, s_mm, s_ss, s_std, amplitude, s_snr))


def decimaldegree_ddmmss(s, dir=None):
    # input:  80.99, dir can be north or east (N or E)
    # output: 805835.4, direction is lost if dir is None
    sign = 1 if s >= 0 else -1
    s = abs(s)
    deg = math.floor(s)
    mint = (s - deg) * 60.0
    min = math.floor(mint)
    sec = (mint - min) * 60.0

    if dir is None:
        r = '%02d%02d%04.1f' % (deg, min, sec)
    else:
        if dir == 'N':
            d = 'N' if sign > 0 else 'S'
            r = '%02d%02d%04.1f%s' % (deg, min, sec, d)
        elif dir == 'E':
            d = 'E' if sign > 0 else 'W'
            r = '%02d%02d%04.1f%s' % (deg, min, sec, d)
        else:
            raise "unknown direction specifier"

    return r


def write_stations_dat(stations, filename):
    # set up stations
    with open(filename, 'w') as sfd:
        for s in stations:
            lon = decimaldegree_ddmmss(s[1], 'E')
            lat = decimaldegree_ddmmss(s[2], 'N')

            llon = ""
            k = 0
            while k < lon.find('.') - 1:
                if lon[k] == '0':
                    llon += " "
                else:
                    break
                k += 1

            llon += lon[k:]
            lon = llon

            llat = ""
            k = 0
            while k < lat.find('.') - 1:
                if lat[k] == '0':
                    llat += " "
                else:
                    break
                k += 1

            llat += lat[k:]
            lat = llat

            sfd.write("{0:<5s} {1:>9}{2:>10} {3:6.1f}\n".format(s[0], lat, lon, s[3]))


def read_hyposat_out(filename):
    event_solution_stats = dict()
    event_solution_stats['magnitude'] = ''
    with open(filename, 'r') as f:
        lines = f.readlines()
        solution = False
        next_line = False
        for line in lines:
            if next_line:
                source_time_UTC = obspy.UTCDateTime(line[:23])
                lat = float(line[23:].split()[0])
                lon = float(line[23:].split()[1])
                z = float(line[23:].split()[2])
                vpvs = float(line[23:].split()[3])
                rms = float(line[23:].split()[-1])
                next_line = False

            if line[:9] == 'Magnitude':
                magnitude = line[11:-1]
                event_solution_stats['magnitude'] = magnitude
            elif line[:30] == 'T0                         LAT':
                next_line = True
                solution = True

    if solution:
        event_solution_stats['source_time_UTC'] = source_time_UTC
        event_solution_stats['lat'] = lat
        event_solution_stats['lon'] = lon
        event_solution_stats['z'] = z
        event_solution_stats['vpvs'] = vpvs
        event_solution_stats['rms'] = rms
        return event_solution_stats
    else:
        return None


