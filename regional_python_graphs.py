import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
import sys
import dropbox_client

putOnDropbox = True

def printRegionalVolume(data, ax, col, ymin, ymax, name):
	regional = data[1:,col]
	regional = np.array([i.lstrip() for i in regional]).astype(float)
	padded = np.pad(regional, (177-171, 171+177*15 - regional.shape[0]), 'constant', constant_values=(np.nan,))
	matrix = padded.reshape((16,177))
	matrix = matrix/1000.0
	dates = np.arange(1,178)
	
	ax.plot(dates, matrix[-16,:], label='2010/11', color=(0.65,0.65,0.65));
	ax.plot(dates, matrix[-15,:], label='2011/12', color=(0.44,0.19,0.63));
	ax.plot(dates, matrix[-14,:], label='2012/13', color=(0.0,0.13,0.38));
	ax.plot(dates, matrix[-13,:], label='2013/14', color=(0,0.44,0.75));
	ax.plot(dates, matrix[-12,:], label='2014/15', color=(0.0,0.69,0.94));
	ax.plot(dates, matrix[-11,:], label='2015/16', color=(0,0.69,0.31));
	ax.plot(dates, matrix[-10,:], label='2016/17', color=(0.57,0.82,0.31));
	ax.plot(dates, matrix[-9,:], label='2017/18', color=(1.0,0.75,0));
	ax.plot(dates, matrix[-8,:], label='2018/19', color=(0.9,0.4,0.05));
	ax.plot(dates, matrix[-7,:], label='2019/20', color=(1.0,0.5,0.5));
	ax.plot(dates, matrix[-6,:], label='2020/21', color=(0.58,0.54,0.33));
	ax.plot(dates, matrix[-5,:], label='2021/22', color=(0.4,0,0.2));
	ax.plot(dates, matrix[-4,:], label='2022/23', color=(0.6,0.6,0.2));
	ax.plot(dates, matrix[-3,:], label='2023/24', color=(0.7,0.2,0.3));
	ax.plot(dates, matrix[-2,:], label='2024/25', color=(0.3,0.2,0.3));
	ax.plot(dates, matrix[-1,:], label='2025/26', color=(1.0,0,0), linewidth=3);
	#ax.set_xlabel("day")
	ax.set_ylabel("Sea ice volume (10$^3\!$ km$^3\!$)")
	ax.set_title(name)
	ax.legend(loc=4, prop={'size': 8})#, bbox_to_anchor=(0.75,1))
	#ax.text(75, .025, 'some text')
	#ax.text(2.5, 2.5, r'$\mu=115,\ \sigma=15$')
	ax.axis([0, 177, ymin, ymax])
	#ax.set_xlim([0, 365])
	#ax.set_ylim([ymin, ymax])
	ax.grid(True);
	#ax.set_xticks(np.arange(0, 1, step=0.2))
	ax.set_xticks(np.arange(3), ['Oct', 'Nov', 'Dec'])  # Set text labels.
	months = ['Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr']
	#print(np.linspace(0,177,7))
	ax.set_xticks([0,14,44,75,106,134,165], ['', '', '', '', '', '', ''])
	ax.xaxis.set_minor_locator(ticker.FixedLocator([7,29,59.5,90.5,120,149.5,171]))
	ax.xaxis.set_minor_formatter(ticker.FixedFormatter(months))
	ax.tick_params(which='minor', length=0)
	#for label in ax.get_xticklabels():
	#	label.set_horizontalalignment('center')
	#for tick in ax.xaxis.get_minor_ticks():
	#	tick.label1.set_horizontalalignment('center')

	
def saveRegionalPlot(col, ymin, ymax, data, name, filename):
	#print('inside saveRegionalPlot', name)
	fig, axs = plt.subplots(figsize=(8, 5))
	printRegionalVolume(data, axs, col, ymin, ymax, name)
	fig.savefig(filename)

def plotRegionalGraphs():
	csvFileName = "cryosat-smos-regional-volume.csv"
	data = np.loadtxt(csvFileName, delimiter=",", dtype=str)

	saveRegionalPlot(4, 0, 1.4, data, "Beaufort Sea CryoSat-SMOS ice volume", "cryosat-smos-volume-beaufort.png")
	saveRegionalPlot(5, 0, 1.4, data, "Chukchi Sea CryoSat-SMOS ice volume", "cryosat-smos-volume-chukchi.png")
	saveRegionalPlot(6, 0, 1.4, data, "East Siberian Sea CryoSat-SMOS ice volume", "cryosat-smos-volume-ess.png")
	saveRegionalPlot(7, 0, 0.7, data, "Laptev Sea CryoSat-SMOS ice volume", "cryosat-smos-volume-laptev.png")
	saveRegionalPlot(8, 0, 1.0, data, "Kara Sea CryoSat-SMOS ice volume", "cryosat-smos-volume-kara.png")
	saveRegionalPlot(9, 0, 0.5, data, "Barents Sea CryoSat-SMOS ice volume", "cryosat-smos-volume-barents.png")
	saveRegionalPlot(10, 0, 1.0, data, "Greenland Sea CryoSat-SMOS ice volume", "cryosat-smos-volume-greenland.png")
	saveRegionalPlot(11, 3, 11, data, "Central Arctic Basin CryoSat-SMOS ice volume", "cryosat-smos-volume-cab.png")
	saveRegionalPlot(12, 0, 1.6, data, "Canadian Arctic Archipelago CryoSat-SMOS ice volume", "cryosat-smos-volume-caa.png")
	saveRegionalPlot(13, 0, 1.0, data, "Baffin Bay CryoSat-SMOS ice volume", "cryosat-smos-volume-baffin.png")
	saveRegionalPlot(14, 0, 1.3, data, "Hudson Bay CryoSat-SMOS ice volume", "cryosat-smos-volume-hudson.png")
	saveRegionalPlot(15, 0, 0.12, data, "Other CryoSat-SMOS ice volume", "cryosat-smos-volume-other.png")
	saveRegionalPlot(16, 0, 21, data, "Total CryoSat-SMOS ice volume", "cryosat-smos-volume-total.png")
	saveRegionalPlot(3, 0, 0.5, data, "Bering Sea CryoSat-SMOS ice volume", "cryosat-smos-volume-bering.png")
	saveRegionalPlot(2, 0, 0.35, data, "Sea of Okhotsk CryoSat-SMOS ice volume", "cryosat-smos-volume-okhotsk.png")
	
	saveRegionalPlot(16, 0, 22, data, "CryoSat-SMOS Arctic sea ice volume", "cryosat-smos-volume.png")

	fig, axs = plt.subplots(7, 2, figsize=(16, 35))
	fig.tight_layout(pad=5.0)

	printRegionalVolume(data, axs[0][0], 4, 0, 1.4, "Beaufort Sea CryoSat-SMOS ice volume")
	printRegionalVolume(data, axs[0][1], 5, 0, 1.4, "Chukchi Sea CryoSat-SMOS ice volume")
	printRegionalVolume(data, axs[1][0], 6, 0, 1.4, "East Siberian Sea CryoSat-SMOS ice volume")
	printRegionalVolume(data, axs[1][1], 7, 0, 0.7, "Laptev Sea CryoSat-SMOS ice volume")
	printRegionalVolume(data, axs[2][0], 8, 0, 1.0, "Kara Sea CryoSat-SMOS ice volume")
	printRegionalVolume(data, axs[2][1], 9, 0, 0.5, "Barents Sea CryoSat-SMOS ice volume")
	printRegionalVolume(data, axs[3][0], 10, 0, 1.0, "Greenland Sea CryoSat-SMOS ice volume")
	printRegionalVolume(data, axs[3][1], 11, 3, 11, "Central Arctic Basin CryoSat-SMOS ice volume")
	printRegionalVolume(data, axs[4][0], 12, 0, 1.6, "Canadian Arctic Archipelago CryoSat-SMOS ice volume")
	printRegionalVolume(data, axs[4][1], 13, 0, 1.0, "Baffin Bay CryoSat-SMOS ice volume")
	printRegionalVolume(data, axs[5][0], 14, 0, 1.3, "Hudson Bay CryoSat-SMOS ice volume")
	printRegionalVolume(data, axs[5][1], 3, 0, 0.5, "Bering Sea CryoSat-SMOS ice volume")
	printRegionalVolume(data, axs[6][0], 2, 0, 0.35, "Sea of Okhotsk CryoSat-SMOS ice volume")
	printRegionalVolume(data, axs[6][1], 16, 0, 21, "Total CryoSat-SMOS ice volume")
	#axs[6][1].axis('off')
	#axs[4][2].axis('off')

	#fig.show()
	#wait = input("Press Enter to continue.")
	dropboxFileName = 'cryosat-smos-regional-volume-7x2.png'
	fig.savefig(dropboxFileName)
	if putOnDropbox:
		dropbox_client.uploadToDropbox([csvFileName, "cryosat-smos-volume-total.png", 'cryosat-smos-thickness-latest.png', 'cryosat-smos-thickness-anomaly-latest.png'])
		#uploadToDropbox(csvFileName)		
		#uploadToDropbox(dropboxFileName)


print('__name__: ',__name__)
if __name__ == "__main__":
	plotRegionalGraphs()