import numpy as np
from struct import unpack
import datetime
import pdb

# define instrument specifications
# FOV: field of view
# pol: polarization at nadir view

instruments = {'amsua':
		  {'name':'AMSU-A',
		  'FOV':3.33,
		  'footprint':[48.,48.,48.,48.,48.,48.,48.,48.,48.,48.,48.,48.,48.,48.,48.,48.],
		  'freq_str':[],
		  'pol':['V','V','V','V','H','H','V','H','H','H','H','H','H','H']},
		'amsub':
		  {'name':'AMSU-B',
		  'FOV':1.1,
		  'footprint':[16.,16.,16.,16.,16.],
		  'freq_str':[],
		  'pol':['V','V','H','H','V']},
		'mhs':
		  {'name':'MHS',
		  'FOV':1.1,
		  'footprint':[16.,16.,16.,16.,16.],
		  'freq_str':['89.0','157.0','183.311$\pm$1.0','183.311$\pm$3.0','190.311'],
		  'pol':['V','V','H','H','V']}}
spacecraft = {'noaa18':
		{'name':'NOAA-18',
		 'height':833.}}

class swath(object):
  def __init__(self):
    pass
  
  def read_swath(self,fname,sensor='amsub',satellite='noaa18'):
    p1 = 1.1910658e-5
    p2 = 1.438833

    self.sensor = sensor
    self.name = instruments[sensor]['name']
    self.fov = instruments[sensor]['FOV']
    self.footprint = instruments[sensor]['footprint']
    self.pol = instruments[sensor]['pol']
    self.height = spacecraft[satellite]['height']
    fid = open(fname, 'rb')
    # read header 
    head = Head()
    head.read_1bhead(fid)
    self.nrec = head.last_scan_record
    nrec = self.nrec
    self.sat_id = head.sat_id
    self.sensor = head.sensor
    sensor = self.sensor
    self.nchan = head.nchan
    nchan = self.nchan
    self.npos = head.npos
    npos = self.npos
    self.tcount = head.tcount
    tcount = self.tcount
    self.start_utc_of_day = head.start_utc_of_day
    self.end_utc_of_day = head.end_utc_of_day
    self.tr_central_wave_number = head.tr_central_wave_number
#        self.hour = (self.start_utc_of_day + 8000 * arange(nrec)) / (1000.*3600.)
    self.scanline_number = np.zeros([nrec],dtype=np.int64)
    self.scanline_year = np.zeros([nrec],dtype=np.int64)
    self.scanline_day_of_year = np.zeros([nrec],dtype=np.int64)
    self.scanline_utc_time_of_day = np.zeros([nrec],dtype=np.int64)
    self.date_time = np.empty([nrec], dtype=object)
    self.mode = np.zeros([nrec],dtype=np.int64)
    self.altitude = np.zeros([nrec])
    self.sza = np.zeros([nrec, npos])
    self.vza = np.zeros([nrec, npos])
    self.lat = np.zeros([nrec, npos])
    self.lon = np.zeros([nrec, npos])
    self.counta = np.zeros([nrec, npos, nchan],dtype=np.int64)
    self.countc = np.zeros([nrec, tcount, nchan],dtype=np.int64)
    self.countw = np.zeros([nrec, tcount, nchan],dtype=np.int64)
    if (sensor == 'A'):
        self.countwd = np.zeros([nrec,17],dtype=np.int64)
        self.countwa = np.zeros([nrec,3],dtype=np.int64)
    self.rada = np.zeros([nrec, npos, nchan])
    self.radc = np.zeros([nrec, tcount, nchan])
    self.radw = np.zeros([nrec, tcount, nchan])
    self.ta = np.zeros([nrec, npos, nchan])
    self.tc = np.zeros([nrec, tcount, nchan])
    self.tw = np.zeros([nrec, tcount, nchan])
    if (sensor == 'A'):
        self.twd = np.zeros([nrec,17])
        self.twa = np.zeros([nrec,3])
    self.clcoef = np.zeros([nrec, 3, nchan])

    if (sensor == 'A'):
        dcoef = np.zeros([4,17])
        acoef = np.zeros([2,3])
        dcoef[0,:] = head.warm_digital_coeffs[0,:] * 1.e-4
        dcoef[1,:] = head.warm_digital_coeffs[1,:] * 1.e-9
        dcoef[2,:] = head.warm_digital_coeffs[2,:] * 1.e-16
        dcoef[3,:] = head.warm_digital_coeffs[3,:] * 1.e-20
        acoef = head.warm_analog_coeffs * 1.e-3

    def convert_to_h_m_s_m(milliseconds):
        """Return the tuple of hours, minutes, seconds, and milliseconds."""
    
        seconds, milliseconds = divmod(milliseconds,1000)
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
#        days, hours = divmod(hours, 24)
    
        return hours, minutes, seconds, milliseconds
    # read data
    for irec in range(nrec):
        scanline = Scanline()
        if (sensor == 'A'):
            scanline.read_amsua_1bdata(fid,irec)
        elif (sensor == 'B') | ( sensor == 'S'):
            scanline.read_amsubmhs_1bdata(fid,irec)
        self.scanline_number[irec] = scanline.scanline_number
        self.scanline_year[irec] = scanline.scanline_year
        self.scanline_day_of_year[irec] = scanline.scanline_day_of_year
        self.scanline_utc_time_of_day[irec] = scanline.scanline_utc_time_of_day
        hours, minutes, seconds, milliseconds = convert_to_h_m_s_m(self.scanline_utc_time_of_day[irec])
        dt = datetime.date(self.scanline_year[0],1,1)+datetime.timedelta(self.scanline_day_of_year[0])
        self.date_time[irec] = datetime.datetime(self.scanline_year[irec],dt.month,dt.day,hours,minutes,seconds,milliseconds)
        self.mode[irec] = scanline.mode >> 15
        self.altitude[irec] = scanline.altitude *1.e-1
        self.sza[irec,:] = scanline.sza * 1.e-2
        self.vza[irec,:] = scanline.vza * 1.e-2
        self.lat[irec,:] = scanline.lat * 1.e-4
        self.lon[irec,:] = scanline.lon * 1.e-4
        self.counta[irec,:,:] = scanline.counta
        self.countc[irec,:,:] = scanline.countc
        self.countw[irec,:,:] = scanline.countw
        if (sensor == 'A'):
            self.countwd[irec,:] = scanline.countwd
            self.countwa[irec,:] = scanline.countwa
            self.clcoef[irec,0,:] = scanline.pri_cal_coeffs[0,:] * 1.e-19
            self.clcoef[irec,1,:] = scanline.pri_cal_coeffs[1,:] * 1.e-13
            self.clcoef[irec,2,:] = scanline.pri_cal_coeffs[2,:] * 1.e-9
        elif (sensor == 'B') | ( sensor == 'S'):
            self.clcoef[irec,0,:] = scanline.pri_cal_coeffs[0,:] * 1.e-16
            self.clcoef[irec,1,:] = scanline.pri_cal_coeffs[1,:] * 1.e-10
            self.clcoef[irec,2,:] = scanline.pri_cal_coeffs[2,:] * 1.e-6
        for ichan in range(nchan):
            self.rada[irec,:,ichan] = (self.clcoef[irec,2,ichan] + 
                self.clcoef[irec,1,ichan] * self.counta[irec,:,ichan] + 
                self.clcoef[irec,0,ichan] * self.counta[irec,:,ichan]** 2)
            self.radc[irec,:,ichan] = (self.clcoef[irec,2,ichan] + 
                self.clcoef[irec,1,ichan] * self.countc[irec,:,ichan] + 
                self.clcoef[irec,0,ichan] * self.countc[irec,:,ichan]** 2)
            self.radw[irec,:,ichan] = (self.clcoef[irec,2,ichan] + 
                self.clcoef[irec,1,ichan] * self.countw[irec,:,ichan] + 
                self.clcoef[irec,0,ichan] * self.countw[irec,:,ichan]** 2)
        if (sensor == 'A'):
                self.twd[irec,:] = (dcoef[0,:] + 
                    dcoef[1,:] * self.countwd[irec,:] + 
                    dcoef[2,:] * self.countwd[irec,:]**2 + 
                    dcoef[3,:] * self.countwd[irec,:]**3)
                self.twa[irec,:] = (acoef[0,:] + 
                    acoef[1,:] * self.countwa[irec,:])

    if (sensor == 'A'): self.twa = self.twa * 2.e-2
    
    fid.close()
    
    for ichan in range(nchan):
        pl = head.tr_central_wave_number[ichan] * 1.e-6
        c1v3 = p1 * pl ** 3
        c2v  = p2 * pl
        tmp  = np.log(c1v3 / self.rada[:,:,ichan] + 1)
        self.ta[:,:,ichan] = c2v / tmp
        tmp  = np.log(c1v3 / self.radc[:,:,ichan] + 1)
        self.tc[:,:,ichan] = c2v / tmp
        tmp  = np.log(c1v3 / self.radw[:,:,ichan] + 1)
        self.tw[:,:,ichan] = c2v / tmp

    return

    def __del__(self):
        pass

class Head:
    def __init__(self):
        pass
    
    def read_1bhead(self, fid):
        fid.seek(534,0)      # if archive header record (512) + 22
#        fid.seek(22, 0)                                         # 22
        self.data_set_name = unpack('>42c',fid.read(42)) # get data set name
        self.sensor = self.data_set_name[6]
        # set sensor specific data
        if (self.sensor == 'A'):
            self.nchan = 15
            self.npos = 30
            self.tcount = 2
            self.warm_digital_coeffs = np.zeros([4,17],dtype=np.int64)
            self.warm_analog_coeffs = np.zeros([2,3],dtype=np.int64)
        elif (self.sensor == 'B') | (self.sensor == 'S'):
            self.nchan = 5
            self.npos = 90
            self.tcount = 4

        self.tr_central_wave_number = np.zeros(self.nchan,dtype=np.int64)

        fid.seek(8,1)
        self.sat_id = unpack('>H',fid.read(2))[0]                     # 74
        self.ins_id = unpack('>2B', fid.read(2))[0]                   # 76
        
        fid.seek(8, 1)                                                # 84
        self.start_year = unpack('>H', fid.read(2))[0]                # 86
        self.start_day_of_year = unpack('>H', fid.read(2))[0]         # 88
        self.start_utc_of_day = unpack('>I', fid.read(4))[0]          # 92

        fid.seek(4, 1)                                                # 96
        self.end_year = unpack('>H', fid.read(2))[0]                  # 98
        self.end_day_of_year = unpack('>H', fid.read(2))[0]           # 100
        self.end_utc_of_day = unpack('>I', fid.read(4))[0]            # 104
        
        if (self.sensor == 'A'):
            fid.seek(40, 1)                                           # 144
            self.last_scan_record = unpack('>H', fid.read(2))[0]      # 146

            fid.seek(542, 1)
            for i in range(self.nchan):                                      # 688
                self.tr_central_wave_number[i] = unpack('>I', fid.read(4))[0]# 868
                fid.seek(8,1)
                
            fid.seek(668, 1)
            temp = np.array(unpack('>40i', fid.read(4*10*4))) # 1536
            self.warm_digital_coeffs[:,0:10] = temp.reshape([10,4]).transpose() # 1696

            fid.seek(36, 1)                                          # 1732
            temp = np.array(unpack('>4i', fid.read(2*2*4)))
            self.warm_analog_coeffs[:,0:2] = temp.reshape([2,2]).transpose()  # 1748

            fid.seek(364, 1)                                         # 2112
            self.warm_digital_coeffs[:,16] = unpack('>i', fid.read(4)) # 2128
            temp = np.array(unpack('>24i', fid.read(4*6*4)))
            self.warm_digital_coeffs[:,10:16] = temp.reshape([6,4]).transpose() # 2224

            fid.seek(28, 1)                                          # 2252
            self.warm_analog_coeffs[:,2] = unpack('>2i', fid.read(8))      # 2260

        elif (self.sensor == 'B'):
            fid.seek(28, 1)
            self.last_scan_record = unpack('>H', fid.read(2))[0]   # 134
        
            fid.seek(190, 1)                                         # 324
            for i in range(self.nchan):
                self.tr_central_wave_number[i] = unpack('>I', fid.read(4))[0]
                fid.seek(8,1)   # 384
                
        elif (self.sensor == 'S'):
            fid.seek(28,1)
            self.last_scan_record = unpack('>H', fid.read(2))[0]   # 134

            fid.seek(282, 1)                                         # 416
            for i in range(self.nchan):
                self.tr_central_wave_number[i] = unpack('>I', fid.read(4))[0]
                fid.seek(8,1)   # 476
class Scanline:
    def __init__(self):
        pass
    
    def read_amsua_1bdata(self,fid,irec):
        L_HBLOCK = 2560
        L_SCANLINE = 2560
        nchan = 15
        npos = 30
        
        self.counta = np.zeros([npos, nchan],dtype=np.int64)                       # antenna
        self.countc = np.zeros([2, nchan],dtype=np.int64)                          # cold target
        self.countw = np.zeros([2, nchan],dtype=np.int64)                          # warm target
        self.countwd = np.zeros(17,dtype=np.int64)                              # warm load digital
                    # 1-4, a1-1; 5, a1-1 center;
                    # 6-9, a1-2; 10, a1-2 center;
                    # 11-16, a2; 17, a2 center;
        self.countwa = np.zeros(3,dtype=np.int64);                           # warm load analog 
                    # 1, a1-1; 2, a1-2; 3, a2;
        
        fid.seek(L_HBLOCK + (irec) * L_SCANLINE + 512, 0)       # 0
#        fid.seek(12, 1)                                          # 12
        self.scanline_number = unpack('>H', fid.read(2))[0]
        self.scanline_year = unpack('>H', fid.read(2))[0]
        self.scanline_day_of_year = unpack('>H', fid.read(2))[0]
        fid.seek(2,1)
        self.scanline_utc_time_of_day = unpack('>I', fid.read(4))[0] # Scan Line UTC Time of Day in milliseconds
        self.mode = unpack('>H', fid.read(2))[0]                    # 14
        fid.seek(66, 1)                                          # 60
        self.pri_cal_coeffs = np.array(unpack('>45i', fid.read(3*nchan*4))).reshape([nchan,3]).transpose() # 120
        
        fid.seek(210, 1)                                          # 210
        self.altitude = unpack('>H', fid.read(2))[0]
        temp = np.array(unpack('>90h', fid.read(npos*3*2)))                         # 752
        self.sza = temp[::3]
        self.vza = temp[1::3]
        
        temp = np.array(unpack('>60i', fid.read(npos*2*4)))                         # 1472
        self.lat = temp[::2]
        self.lon = temp[1::2]

        fid.seek(12, 1)                                           # 1480
        temp = np.array(unpack('>510H', fid.read(1020)))                           # 2560
        for ichan in range(2,nchan):
            self.counta[:,ichan] = temp[ichan+2::17]
        
        temp = np.array(unpack('>30H', fid.read(60)))                            # 2616
        self.countc[0,2:15] = temp[4:17]
        self.countc[1,2:15] = temp[17::]
        
        temp = np.array(unpack('>46H', fid.read(92)))                            # 2616
        self.countwd[0:10] = temp[35:45]
        
        temp = np.array(unpack('>30H', fid.read(60)))                            # 2664
        self.countw[0,2:15] = temp[4:17]
        self.countw[1,2:15] = temp[17::]

        fid.seek(16, 1)                                          # 2152
        temp = unpack('>28B', fid.read(28))                             # 2180
        self.countwa[0:2] = temp[4:6]

        fid.seek(12, 1)                                          # 2192
        temp = unpack('>120H', fid.read(240))                           # 2432
        self.counta[:,0] = temp[0+2::4]
        self.counta[:,1] = temp[1+2::4]

        temp = unpack('>6H', fid.read(12))                             # 2444
        self.countc[0,0:2] = temp[2:4]
        self.countc[1,0:2] = temp[4::]

        temp = unpack('>20H', fid.read(40))                            # 2484
        self.countwd[16] = temp[12]
        self.countwd[10:16] = temp[13:19]

        temp = unpack('>6H', fid.read(12))                             # 2496
        self.countw[0,0:2] = temp[2:4]
        self.countw[1,0:2] = temp[4:6]

        fid.seek(19, 1)                                          # 2515
        self.countwa[2] = unpack('>B', fid.read(1))[0]
        
    def read_amsubmhs_1bdata(self,fid,irec):
        L_HBLOCK = 3072
        L_SCANLINE = 3072
        nchan = 5
        npos = 90
        
        self.counta = np.zeros([npos, nchan],dtype=np.int64)                       # antenna
        self.countc = np.zeros([4, nchan],dtype=np.int64)                          # cold target
        self.countw = np.zeros([4, nchan],dtype=np.int64)                          # warm target
        self.lat = np.zeros(npos,dtype=np.int64)
        self.lon = np.zeros(npos,dtype=np.int64)
        
        fid.seek(L_HBLOCK + (irec) * L_SCANLINE + 512, 0)       # 0
#        fid.seek(12, 1)                                          # 12
        self.scanline_number = unpack('>H', fid.read(2))[0]
        self.scanline_year = unpack('>H', fid.read(2))[0]
        self.scanline_day_of_year = unpack('>H', fid.read(2))[0]
        fid.seek(2,1)
        self.scanline_utc_time_of_day = unpack('>I', fid.read(4))[0]
        self.mode = unpack('>H', fid.read(2))[0]                    # 14
        fid.seek(46, 1)                                          # 60
        self.pri_cal_coeffs = np.array(unpack('>15i', fid.read(3*nchan*4))).reshape([nchan,3]).transpose() # 120
        
        fid.seek(90, 1)                                          # 210
        self.altitude = unpack('>H', fid.read(2))[0]
        temp = np.array(unpack('>270h', fid.read(npos*3*2)))                         # 752
        self.sza = temp[::3]
        self.vza = temp[1::3]
        
        temp = np.array(unpack('>180i', fid.read(npos*2*4)))                         # 1472
        self.lat = temp[::2]
        self.lon = temp[1::2]

        fid.seek(8, 1)                                           # 1480
        temp = np.array(unpack('>540H', fid.read(1080)))                           # 2560
        for ichan in range(nchan):
            self.counta[:,ichan] = temp[ichan+1::6]
        
        
        fid.seek(8, 1)                                           # 2568
        temp = np.array(unpack('>24H', fid.read(48)))                            # 2616
        self.countc[0,:] = temp[1:6]
        self.countc[1,:] = temp[7:12]
        self.countc[2,:] = temp[13:18]
        self.countc[3,:] = temp[19::]
        
        temp = np.array(unpack('>24H', fid.read(48)))                            # 2664
        self.countw[0,:] = temp[1:6]
        self.countw[1,:] = temp[7:12]
        self.countw[2,:] = temp[13:18]
        self.countw[3,:] = temp[19::]


    def __del__(self):
        pass    
  
