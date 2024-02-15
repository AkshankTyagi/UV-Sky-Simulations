pro read_castelli, wave, data, itemper
star_dir = "/Users/jayanth/user/data/castelli/ckp00/"
star_dir = 'C:\Users\jmurt\Dropbox\user\data\castelli\ckp00\'

;Define temperatures and g
	temper=[50000,49000,48000,47000,46000,45000,44000,43000,42000,41000,40000]
	gindex=[11,    11,   11,   11,  11,   10,   10,   10, 10,      10,   10]
	temper=[temper,39000,38000,37000,36000,35000,34000,33000,32000,31000,30000]
	gindex=[gindex, 9,  9,    9,   9,     9,    9,    9,    9,    9,    9]
	temper=[temper,29000,28000,27000,26000,25000,24000,23000,22000,21000,20000]
	gindex=[gindex, 9,    9,    9,   9,     9,    9,    9,    9,    9,    9]
	temper=[temper,19000,18000,17000,16000,15000,14000,13000,12750,12500]
	gindex=[gindex, 9,    9,    9,     9,     9,    9,    9,    9,    9]
	temper=[temper,12250,12000,11750,11500,11250,11000,10750,10500,10250,10000]
	gindex=[gindex, 9,  9,    9,   9,     9,    9,    9,    9,    9,    9]
	temper=[temper,9750,9500,9250,9000,8750,8500,8250,8000,7750,7500,7250,7000]
	gindex=[gindex, 9,  9,    9,   9,     9,    9,  9,   9,  9,    9,  9,   9]
	temper=[temper,6750,6500,6250,6000,5750,5500,5250,5000,4750,4500,4250,4000]
	gindex=[gindex, 9,  9,    10,   10,  10,  10,  10,  10,   10,  10, 10, 10]
	temper=[temper,3750,3500]
	gindex=[gindex, 10,  10]
	filename = star_dir + "ckp00_"+strtrim(string(temper(itemper)),1)+".fits"
	im=mrdfits(filename,1,hdr,/silent)
	ck_wave=im.wavelength
	ck_data=im.(gindex(itemper))
	if (n_elements(wave) eq 0)then begin
		wave = ck_wave
		data = ck_data
	endif else data = interpol(ck_data, ck_wave, wave)
	stop
end