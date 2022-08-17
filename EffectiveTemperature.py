import numpy as np
import os
import matplotlib.pyplot as plt
from datetime import timezone,date,timedelta,datetime
from lmfit import Model
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.stats import gaussian_kde

def get_season(now):
	Y = 2000
	seasons = [('winter', (date(Y,  1,  1),  date(Y,  3, 20))), ('spring', (date(Y,  3, 21),  date(Y,  6, 20))), ('summer', (date(Y,  6, 21),  date(Y,  9, 22))), ('autumn', (date(Y,  9, 23),  date(Y, 12, 20))), ('winter', (date(Y, 12, 21),  date(Y, 12, 31)))]
	if isinstance(now, datetime): now = now.date()
	now = now.replace(year=Y)
	return next(season for season, (start, end) in seasons if start <= now <= end)

def do_unique(arr_day,arr_temp,arr_press,eff_Temp=np.array([])):
	day = arr_day[0]
	aux_indx = np.where(arr_day==day)[0]

	mean_eff_temp = np.array([ get_eff_temp(arr_temp[i],arr_press[i],E_EH=energy) for i in aux_indx]).mean()

	aux_day = np.delete(arr_day,aux_indx)
	aux_temp = np.delete(arr_temp,aux_indx)
	aux_press = np.delete(arr_press,aux_indx)

	aux_effTemp = eff_Temp
	aux_effTemp = np.append(aux_effTemp,mean_eff_temp)

	if aux_day.size==0: return aux_effTemp

	return do_unique(aux_day,aux_temp,aux_press,aux_effTemp)

def gaussian(x, amp, cen, wid):
    """1-d gaussian: gaussian(x, amp, cen, wid)"""
    return (amp / (np.sqrt(2*np.pi) * wid)) * np.exp(-(x-cen)**2 / (2*wid**2))

def lbd(p):
	if p=="pi": return 180.
	return 160.

def l(p):
	#
	return 1./(1./120. - 1./lbd(p))

def eps(p,e_pi,e_K):
	if p=="pi": return e_pi
	return e_K

def A(p,r_K):
	if p=="pi": return 1.
	return r_K*0.38

def B(p,b_pi,b_K):
	if p=="pi": return b_pi
	return b_K

def K(x,p):
	aux1 = x*((1-x/l(p))**2)
	aux2 = (1 - np.exp(-x/l(p)))*l(p)
	return aux1/aux2

def weight(x,p,E_EH,gm,b_pi, b_K, r_K, e_pi, e_K): #p is the particle "pi" or "K"
	aux1 = (1 - x/l(p))**2
	aux2 = np.exp(-x/lbd(p))*A(p,r_K)
	aux3 = gm + (gm + 1.)*B(p,b_pi,b_K)*K(x,p)*((E_EH/eps(p,e_pi,e_K))**2)
	return (aux1*aux2)/aux3

def get_eff_temp(atemp, apress, E_EH, gm=1.7, b_pi=1.46, b_K=1.74, r_K=0.149, e_pi=114., e_K=851.):
	x = apress*1.019 #g/cm2
	dx = np.array([x[i]-x[i-1] for i in range(1,len(x))])
	# dx = np.insert(dx,0,dx.mean())
	x = x[:-1]
	W = weight(x,"pi", E_EH, gm, b_pi, b_K, r_K, e_pi, e_K) + weight(x,"K", E_EH, gm, b_pi, b_K, r_K, e_pi, e_K) #array-like
	n = dx*atemp[:-1]*W
	d = dx*W
	return np.sum(n)/np.sum(d)

def Igra_data(file):
	first_day = (date(2011,12,24)-date(1970,1,1))/timedelta(seconds=1)
	last_day = (date(2017,8,30)-date(1970,1,1))/timedelta(seconds=1)
	on_date = False
	bad_values = [-8888,-9999]
	temp_tmp = np.array([])
	temp_press = np.array([])
	previous_day = 0

	dates = []
	all_temperatures = []
	all_pressures = []

	for line in file:
		if "#" in line:
			ts = (date(int(line[13:17]),int(line[18:20]),int(line[21:23]))-date(1970,1,1))/timedelta(seconds=1)

			if not(on_date):
				if ts<first_day: continue
				on_date = True

			if previous_day == ts:
				same_day = True
				continue

			same_day = False

			if len(temp_tmp)>13:
				all_temperatures.append(temp_tmp[::-1])
				all_pressures.append(temp_press[::-1])
				dates.append((previous_day - previous_day%(3600*24))/(3600*24))

			previous_day = ts
			temp_tmp = np.array([])
			temp_press = np.array([])
			if ts>last_day: break

			continue

		if not(on_date): continue

		if line[1] == "3": continue

		pressure = int(line[9:15])
		temperature = int(line[22:27])
		if pressure in bad_values or temperature in bad_values: continue

		pressure /= 100. #hPa
		if not(pressure in [1,2,3,5,7,10,20,30,50,70,100,125,150,175,200,250,300,350,400,450,500,550,600,650,700,750,775,800,825,850,875,900,925,950,975,1000]): continue
		temperature = temperature/10. + 273.15 #Kelvin

		if same_day:
			indx_temp = np.where(temp_press==pressure)[0]
			if indx_temp.size!=0:
				temp_tmp[indx_temp] = (temp_tmp[indx_temp] + temperature)/2
			else:
				insert_indx = np.where(temp_press>pressure)[0]
				if insert_indx.size == 0: insert_indx = 0 #is bigger than the rest, first position
				else: insert_indx = insert_indx[-1]+1 #Descending order
				temp_tmp = np.insert(temp_tmp,insert_indx,temperature)
				temp_press = np.insert(temp_press,insert_indx,pressure)
			continue

		temp_tmp = np.append(temp_tmp,temperature)
		temp_press = np.append(temp_press,pressure)

	return np.array(dates),all_temperatures,all_pressures

def ECMWF_data(filename, E_EH,avail_press=None, gm=1.7, b_pi=1.46, b_K=1.74, r_K=0.149, e_pi=114., e_K=851.,save=True):
	eff_temp = np.array([])
	days = np.array([])
	temps = []
	pressures = []

	all_data = np.genfromtxt("../temperature/"+filename,skip_header=2,skip_footer=10)
	dates = all_data[:,0].astype(int)
	hours = all_data[:,1].astype(int)
	press = all_data[:,2]
	temp = all_data[:,3] #there are 4 values, first one is the closest one

	press *= 100.
	R = 8.31432
	M0 = 28.9644
	g0 = 9.80665
	ientry = 0
	press /= 100

	init = dates[0]
	yr = int((init - init%10000)/10000)
	aux = init - yr*10000
	day = aux%100
	mth = int((aux-day)/100)

	current_date = (date(yr,mth,day)-date(1970,1,1))/timedelta(seconds=1)

	end = dates[-1]
	yr = int((end - end%10000)/10000)
	aux = end - yr*10000
	day = aux%100
	mth = int((aux-day)/100)

	last_day = (date(yr,mth,day)-date(1970,1,1))/timedelta(seconds=1)

	while current_date<=last_day:
		aux_temp = np.array([])
		aux = int(datetime.utcfromtimestamp(current_date).strftime("%Y%m%d"))
		if avail_press != None:
			aux_day = (current_date - current_date%(3600*24))/(3600*24)
			indx_day = np.where(avail_press[0]==aux_day)[0]
			if indx_day.size == 0:
				current_date += 3600*24
				continue
			indx_day = indx_day[0]

		for hr in [0,600,1200,1800]:
			indx_1 = np.where(dates==aux)
			indx_2 = np.where(hours==hr)
			indx = np.intersect1d(indx_1,indx_2)

			use_temp = temp[indx]
			use_press = press[indx]
			if avail_press != None:
				indx_3 = np.intersect1d(use_press,avail_press[1][indx_day],return_indices=True)[1]
				use_temp = use_temp[indx_3]
				use_press = use_press[indx_3]

			if hr == 0:
				temps.append(use_temp)
				pressures.append(use_press)
			aux_temp = np.append(aux_temp,get_eff_temp(use_temp,use_press,E_EH, gm, b_pi, b_K, r_K, e_pi, e_K))

		eff_temp = np.append(eff_temp,aux_temp.mean())
		days = np.append(days,(current_date - current_date%(3600*24))/(3600*24))

		current_date += 3600*24

	return days,eff_temp,temps,pressures

def Int_EffTemp(press,temps,E_EH,gm=1.7,b_pi=1.46,b_K=1.74,r_K=0.149,e_pi=114.,e_K=851.,kind="linear"):
	effTemps_arr = np.array([])
	for i in range(len(press)):
		temp_func = interp1d(x=press[i]*1.096,y=temps[i],kind=kind,fill_value="extrapolate")

		tot_weight = lambda X: weight(X,"pi", E_EH, gm, b_pi, b_K, r_K, e_pi, e_K) + weight(X,"K", E_EH, gm, b_pi, b_K, r_K, e_pi, e_K)
		to_integrate = lambda X: temp_func(X)*(tot_weight(X))
		integration1 = quad(to_integrate,min(press[i])*1.096,max(press[i]*1.096))[0]
		integration2 = quad(tot_weight,min(press[i])*1.096,max(press[i])*1.096)[0]

		effTemps_arr = np.append(effTemps_arr,integration1/integration2)
	return effTemps_arr

nEH = 3 #1,2,3
if nEH == 1: energy = 37.
if nEH == 2: energy = 41.
if nEH == 3: energy = 143.
ecmwf_file = "dayabay2_temps.dat"
nEH = str(nEH)

f = open("../temperature/CHM00059316-data.txt")
Sht_days,Sht_temps,Sht_press = Igra_data(f)
f.close()
f = open("../temperature/CHM00059280-data.txt")
Qcy_days,Qcy_temps,Qcy_press = Igra_data(f)
f.close()

Sht_effTemps = np.array([get_eff_temp(Sht_temps[i],Sht_press[i], E_EH=energy) for i in range(len(Sht_days))])
Qcy_effTemps = np.array([get_eff_temp(Qcy_temps[i],Qcy_press[i], E_EH=energy) for i in range(len(Qcy_days))])
db_days,db_effTemps,db_temps,db_press = ECMWF_data(ecmwf_file, E_EH=energy)

plt.figure()
plt.plot(db_days,db_effTemps,".r",label="EMCWF")
xtks = np.array([(date(2012,1,1)-date(1970,1,1))/timedelta(seconds=1),
		(date(2013,1,1)-date(1970,1,1))/timedelta(seconds=1),
		(date(2014,1,1)-date(1970,1,1))/timedelta(seconds=1),
		(date(2015,1,1)-date(1970,1,1))/timedelta(seconds=1),
		(date(2016,1,1)-date(1970,1,1))/timedelta(seconds=1),
		(date(2017,1,1)-date(1970,1,1))/timedelta(seconds=1),])
xtks = (xtks - xtks%(3600*24))/(3600*24)
plt.xticks(xtks, ["01 Jan 12", "01 Jan 13", "01 Jan 14", "01 Jan 15", "01 Jan 16", "01 Jan 17"])
plt.ylabel("Temperature [K]",horizontalalignment="right",y=1.0)
plt.show()

dT = Sht_effTemps - db_effTemps

plt.figure()
values,edges,patches = plt.hist(dT,bins=np.arange(-5,5+0.08,0.08),histtype="step")
centers = np.array([ (edges[i+1]+edges[i])/2 for i in range(len(edges)-1)])

gmodel = Model(gaussian)
result = gmodel.fit(values, x=centers, cen=1, amp=max(values), wid=0.4)

sig_cen = np.sqrt(result.covar[1][1])
sig_wid = np.sqrt(result.covar[2][2])

plt.plot(centers,values,".r")
plt.plot(centers,result.best_fit,"r")
plt.xlabel(r"$\Delta$T",horizontalalignment='right',x=1.0)
plt.ylabel("Number of Days",horizontalalignment='right',y=1.0)
ax = plt.gca()
ax.text(0.75,0.95,r"$\mu$ = {0:.2f} $\pm$ {1:.2f}".format(result.best_values["cen"],sig_cen),transform=ax.transAxes)
ax.text(0.75,0.90,r"$\sigma$ = {0:.2f} $\pm$ {1:.2f}".format(result.best_values["wid"],sig_wid),transform=ax.transAxes)
plt.xlim(-2.5,2.5)
plt.show()

# Bootstrap --------------------------------------------

nro_rep = int(1e3)

E_EH = np.random.normal(energy,energy*0.07,nro_rep)
gm = np.random.normal(1.7,0.1,nro_rep)
b_pi = np.random.normal(1.46,0.007,nro_rep)
b_K = np.random.normal(1.74,0.028,nro_rep)
r_K = np.random.normal(0.149,0.06,nro_rep)
e_pi = np.random.normal(114.,3.,nro_rep)
e_K = np.random.normal(851.,14.,nro_rep)

print("Done producing {} parameters".format(nro_rep))

plt.figure()
plt.hist(gm,bins=50,histtype="step")
plt.axvline(x=1.7,ls="--",color="r",label=r"$\gamma$ = 1.7")
plt.ylabel("Counts",horizontalalignment="right",y=1.0)
plt.xlabel(r"$\gamma$",horizontalalignment="right",x=1.0)
plt.legend()
plt.show()

db_days,db_effTemps,db_temps,db_press = ECMWF_data(ecmwf_file, E_EH=energy)
dTs = np.array( [ [ db_effTemps [j] - get_eff_temp(db_temps[j],db_press[j],E_EH=E_EH[i], gm=gm[i], b_pi=b_pi[i], b_K=b_K[i], r_K=r_K[i], e_pi=e_pi[i], e_K=e_K[i]) for i in range(nro_rep) ] for j in range(len(db_days)) ] )

sigmas = np.array([])
for j in range(len(db_days)):
	dTs = np.array([db_effTemps[j] - get_eff_temp(db_temps[j],db_press[j],E_EH=E_EH[i], gm=gm[i], b_pi=b_pi[i], b_K=b_K[i], r_K=r_K[i], e_pi=e_pi[i], e_K=e_K[i]) for i in range(nro_rep)])
	perc  = np.percentile(dTs,[16,50,84])
	wid1 = perc[2]-perc[1]
	wid2 = perc[1]-perc[0]
	mean_wid = (wid1+wid2)/2
	sigmas = np.append(sigmas,mean_wid)
	print("Day {0}/{1} sig = {2} mean = {3}".format(j+1,len(db_days)+1,dTs.std(),perc[1]))

Qcy_press = np.array(Qcy_press)
Qcy_temps = np.array(Qcy_temps)
db_press = np.array(db_press)
db_temps = np.array(db_temps)

t_date = np.array([date.fromtimestamp(int(i)*3600*24) for i in Qcy_days])
ssns = np.array([get_season(i) for i in t_date])
summer = np.where(ssns=="summer")
winter = np.where(ssns=="winter")

plt.figure()

t_date = np.array([date.fromtimestamp(int(i)*3600*24) for i in db_days])
ssns = np.array([get_season(i) for i in t_date])
summer = np.where(ssns=="summer")
winter = np.where(ssns=="winter")
for issn in [["Summer","red"],["Winter","blue"]]:
	indx = np.where(ssns==issn[0].lower())


db_av = np.array([])
db_press_lvl = np.sort(np.unique(np.concatenate(db_press)))
w = weight(db_press_lvl*1.019,"pi", E_EH=energy, gm=1.7, b_pi=1.46, b_K=1.74, r_K=0.149, e_pi=114., e_K=851.) + weight(db_press_lvl,"K", E_EH=energy, gm=1.7, b_pi=1.46, b_K=1.74, r_K=0.149, e_pi=114., e_K=851.)
for j in db_press_lvl:
	temp_sum = 0
	count = 0
	for iday in range(len(db_days)):
		ind = np.where(db_press[iday]==j)[0]
		if ind.size!=0:
			temp_sum += db_temps[iday][ind]
			count+=1
	db_av = np.append(db_av,temp_sum/count)

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ln1 = ax1.plot(db_press_lvl, db_av,color=issn[1],label="Temperature",ls="-")

ax2 = ax1.twinx()
ln2 = ax2.plot(db_press_lvl,w,"--",label="Weight")
ax2.set_ylabel("Weight")

ax1.set_xscale("log")
ax1.set_xlabel("Pressure [hPa]")
ax1.set_ylabel("Temperature [K]")
lns = ln1 + ln2
labs = [l.get_label() for l in lns]
ax1.legend(lns,labs,loc=1)

perc  = np.percentile(sigmas,[16,50,84])

mean_Teff = db_effTemps.mean()
sig = np.sqrt(result.best_values["wid"]**2 + perc[1]**2)
print(mean_Teff,sig)
np.savetxt("../temperature/Teff_EH"+nEH+".dat",np.column_stack((db_days,db_effTemps)),header="{0} {1}".format(mean_Teff,sig),fmt="%d %f")
