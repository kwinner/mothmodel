from csv import reader
import numpy as np
from matplotlib import pyplot as plt

# prior = True
# num_runs = 1000

prior = False
num_samples = 200


N = 30
t = range(0,51,5)

to_plot = 1

postfix = '_posterior_final'
if prior:
	postfix = '_prior'


cum_vals = np.zeros((num_samples,2*N))
with open('latent_variable_cum_vals'+postfix+'.csv','rb') as csvfile:
	csvreader = reader(csvfile)
	for n, row in enumerate(csvreader):
		cum_vals[n,:] = np.array([int(v) for v in row])

times = np.zeros((num_samples,2*N))
with open('latent_variable_times'+postfix+'.csv','rb') as csvfile:
	csvreader = reader(csvfile)
	for n, row in enumerate(csvreader):
		times[n,:] = np.array([float(v) for v in row])

ys = np.zeros((num_samples,len(t)))
with open('latent_variable_ys'+postfix+'.csv','rb') as csvfile:
	csvreader = reader(csvfile)
	for n, row in enumerate(csvreader):
		ys[n,:] = np.array([float(v) for v in row])		

with open('latent_variable_mean_ps_ks'+postfix+'.csv','rb') as csvfile:
	csvreader = reader(csvfile)
	ps = np.array([float(v) for v in csvreader.next()])
	ks = np.array([int(v) for v in csvreader.next()])	

with open('latent_variable_mean_ps_zonn_ks'+postfix+'.csv','rb') as csvfile:
	csvreader = reader(csvfile)
	ps_zonn = np.array([float(v) for v in csvreader.next()])
	ks_zonn = np.array([int(v) for v in csvreader.next()])		

if to_plot == 1:
	num_samples_to_plot = 1
	colors = 'bg'
	# for n in range(num_samples_to_plot):
	n=3
	latent_draw, = plt.step(times[n,:],cum_vals[n,:],'b-',alpha=.5,where='post',label='Posterior sample')
	print type(latent_draw)
	# plt.plot(ks,N*ps,'r-',linewidth=2)
	pdf_zonn, = plt.plot(ks_zonn,N*ps_zonn,'r-',linewidth=2,label='Zonneveld')
	obs_points, = plt.plot(t,ys[0,:],'ko', label='Observed counts')
	plt.xlim((min(t)-1,max(t)+1))
	plt.ylim((0,20))
	plt.xlabel('Days after May 1st')
	plt.ylabel('Abundance')
	plt.legend(handles=[latent_draw,pdf_zonn,obs_points])
	plt.title('Step function of abundance over time, N = %d' % N)
	plt.savefig('1_draw.pdf',format='pdf')
elif to_plot == 2:

	latent_draw, = plt.plot(0,0,'b-',label='Posterior sample')
	for n in range(cum_vals.shape[0]):
		plt.step(times[n,:],cum_vals[n,:],'b-',alpha=.04,where='post')
	
	# pdf, = plt.plot(ks,N*ps,'r-',linewidth=2,label='pdf')
	pdf_zonn, = plt.plot(ks_zonn,N*ps_zonn,'r-',linewidth=2,label='Zonneveld')
	obs_points, = plt.plot(t,ys[0,:],'ko',fillstyle='none',mew=2, ms=15, label='Observed counts')
	plt.legend(handles=[latent_draw,pdf_zonn,obs_points])
	plt.xlim((min(t)-1,max(t)+1))
	plt.ylim((0,20))
	plt.xlabel('Days after May 1st')
	plt.ylabel('Abundance')
	plt.title('Step function of abundance over time, %d draws, N = %d, ' % (num_samples,N))
	plt.savefig('200_draws.pdf',format='pdf')
elif to_plot == 3:

	t_plot = range(-5,max(t)+1)
	lower = np.zeros(len(t_plot))
	upper = np.zeros(len(t_plot))
	for i,time in enumerate(t_plot):
		argsort = [np.argsort(time - times[s,time-times[s,:]<0]) for s in range(num_samples)]
		for a in range(len(argsort)):
			if argsort[a].shape == (0,):
				argsort[a] = np.array([0])

		sorted_vals = np.sort(cum_vals[range(num_samples),[a[0] for a in argsort]])
		lower[i] = sorted_vals[num_samples*.025]
		upper[i] = sorted_vals[num_samples*.975]

	plt.plot(t_plot,lower,'b-')
	plt.plot(t_plot,upper,'c-')
	plt.plot(ks,N*ps,'r-',linewidth=2)
	plt.plot(t,ys[0,:],'ko')
	plt.xlim((min(t)-1,max(t)+1))
	plt.title('N = %d' % N)

plt.show()



