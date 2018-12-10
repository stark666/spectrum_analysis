#coding=utf-8
#AWG芯片测试数据处理脚本

#ACPL=all channel pure loss
#修改中心波长计算方式为每个峰中心波长对应的通带
#添加计算通带插损
#对目前数据及系统进行修改 20181019
#添加中心波长，通道间隔，同通道串扰改为'--' 20181021
#添加注释 20181021
import csv,os,math,re,time
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt

class spectrum():                            #光谱类
	def __init__(self,path,ase_file=None):
		self.path=path                      #芯片光谱文档路径
		self.ase_file=ase_file				#宽谱光源光谱文档路径

	def skiprows(self):                     #跳过光谱仪数据数据刚开始的非数据行
		with open(self.path,'rb') as f:
			for index,line in enumerate(f.readlines()):
				row=line.decode('utf-8').split('\t')
				if row[0][0]=='1':
					skiprows=index
					break
		return skiprows

	def read(self):           #to array     #读光谱数据
		skrow = self.skiprows()
		tmp = np.loadtxt(open(self.path,"rb"),delimiter=",",skiprows=skrow)
		return tmp

	def pure_loss(self):                    #芯片光谱数据与宽谱光源数据相减，得到净芯片光谱
		ase_file=spectrum(self.ase_file)
		ase_spec = ase_file.read()
		pure_loss_all=ase_spec[:,0]    #存波长
		pure_loss_all.shape=(pure_loss_all.shape[0],1)
		awg_spec=self.read()
		pure_loss=awg_spec[:,1]-ase_spec[:,1]
		pure_loss.shape=(pure_loss.shape[0],1)
		pure_loss_all=np.concatenate((pure_loss_all,pure_loss), axis=1)    #拼接波长和损耗
		# print self.path
		tmp=re.match('.*?(\d{1,2}).csv',self.path.lower())
		col=tmp.group(1)
		# print col
		df_pl=pd.DataFrame(pure_loss_all,columns=['WL']+['CH '+str(col)])
		df_pl=df_pl.set_index('WL')
		# self.df_pl=df_pl
		return df_pl      #pure loss,to dataframe

	def find_peak_array(self):                     #返回中心波长波峰及对应插损
		loss=self.pure_loss()
		# loss=self.df_pl
		# print loss
		wavelength=loss.index
		peak_loss=np.max(loss)
		# print peak_loss
		peak_index=np.where(loss==peak_loss)[0][0]
		# peak_index=peak_index[0][0]
		# print wavelength[peak_index]
		peak_wavelength=wavelength[peak_index]
		return np.array([peak_wavelength,round(peak_loss,3)])

	def bandwidth(self,depth=1):                 #返回带宽，默认为1dB带宽，由于光谱仪采点有限添加插值法取数
		df=self.pure_loss()
		loss=df.iloc[:,0]
		peak_wl,peak_loss=self.find_peak_array()
		for index,los in enumerate(loss):
			if los<=peak_loss-depth:
				# print wl,los
				los_ll=los
				wl_ll=index

			else:
				los_lr=los
				wl_lr=index

				break
		wl_ll=loss.index[wl_ll]
		wl_lr=loss.index[wl_lr]
		loss_tmp=loss[df.index>=peak_wl]
		# print loss
		for index,los in enumerate(loss_tmp):
			if los>=peak_loss-depth:
				los_rl=los
				wl_rl=index
			else:
				los_rr=los
				wl_rr=index
				break

		wl_rl=loss_tmp.index[wl_rl]
		wl_rr=loss_tmp.index[wl_rr]
		# print '\n\n\n\n'
		print(wl_ll,wl_lr,los_ll,los_lr)
		print(wl_rl,wl_rr,los_rl,los_rr)
		wl_l=(peak_loss-depth-los_ll)/(los_lr-los_ll)*(wl_lr-wl_ll)+wl_ll
		wl_r=(peak_loss-depth-los_rr)/(los_rl-los_rr)*(wl_rr-wl_rl)+wl_rl
		return round(wl_r-wl_l,3)

	def region_interpolation(self,left,right):                 #插值法取数，left,right为左右波长边界，返回值为left到right范围内的光谱数据
		# print left,right
		df=self.pure_loss()
		for wl in df.index:
			if wl<=left:
				wl_ll=wl
			else:
				wl_lr=wl
				break

		los_ll=df.iloc[:,0][wl_ll]
		los_lr=df.iloc[:,0][wl_lr]

		# print wl_ll,los_ll,wl_lr,los_lr
		df_tmp=df[df.index>=(left+right)/2]
		# print loss
		for wl in df_tmp.index:
			if  wl<=right:
				wl_rl=wl
			else:
				wl_rr=wl
				break
		los_rl=df.iloc[:,0][wl_rl]
		los_rr=df.iloc[:,0][wl_rr]

		los_l=(left-wl_ll)/(wl_lr-wl_ll)*(los_lr-los_ll)+los_ll
		los_r=(right-wl_rl)/(wl_rr-wl_rl)*(los_rr-los_rl)+los_rl

		# return los_l,los_r
		df=df[df.index<=right]
		df=df[df.index>=left]
		# print df
		df=df.reset_index()
		# print df
		tmp=pd.DataFrame([[left,round(los_l,3)],[right,round(los_r,3)]],columns=df.columns)
		df=df.append(tmp)
		df=df.set_index('WL')
		df=df.sort_index()
		# print df
		return df


def ACPL(path_list,ase):          #ACPL=all channel pure loss
	tmp=[]
	# print path_list
	for spec in path_list:
		# spec=os.path.join(path,spec)
		tmp.append(spectrum(spec,ase).pure_loss())
	ACPL=pd.concat(tmp,axis=1)           #ACPL=all channel pure loss
	return ACPL


def find_peak_array(series):                     #返回中心波长波峰及对应插损
	loss=series
	wavelenth=series.index
	peak_loss=np.max(loss)
	peak_index=np.where(loss==peak_loss)
	peak_wavelenth=wavelenth[peak_index].tolist()[0]
	return np.array([peak_wavelenth,round(peak_loss,3)])

def bandwidth(series,peak_wl,peak_loss,depth=1):			#返回带宽，默认为1dB带宽，由于光谱仪采点有限添加插值法取数
		# df=self.pure_loss()
		loss=series
		# peak_wl,peak_loss=self.find_peak_array()
		for index,los in enumerate(loss):
			if los<=peak_loss-depth:
				# print wl,los
				los_ll=los
				wl_ll=index
			else:
				los_lr=los
				wl_lr=index

				break
		wl_ll=loss.index[wl_ll]
		wl_lr=loss.index[wl_lr]
		loss_tmp=loss[series.index>=peak_wl]
		# print loss_tmp
		for index,los in enumerate(loss_tmp):
			if los>=peak_loss-depth:
				los_rl=los
				wl_rl=index
			else:
				los_rr=los
				wl_rr=index
				break

		wl_rl=loss_tmp.index[wl_rl]
		wl_rr=loss_tmp.index[wl_rr]
		# print '\n\n\n\n'
		# print wl_ll,wl_lr,los_ll,los_lr
		# print wl_rl,wl_rr,los_rl,los_rr
		wl_l=(peak_loss-depth-los_ll)/(los_lr-los_ll)*(wl_lr-wl_ll)+wl_ll
		wl_r=(peak_loss-depth-los_rr)/(los_rl-los_rr)*(wl_rr-wl_rl)+wl_rl
		return round(wl_r-wl_l,3)

def region_interpolation(series,left,right):			#插值法取数，left,right为左右波长边界，返回值为left到right范围内的光谱数据
	# print left,right
	# df=self.pure_loss()
	for wl in series.index:

		if wl<=left:
			wl_ll=wl
		else:
			wl_lr=wl
			break

	los_ll=series[wl_ll]
	los_lr=series[wl_lr]

	# print(wl_ll,los_ll,wl_lr,los_lr)
	se_tmp=series[series.index>=(left+right)/2]
	# print('xixi')
	for wl in se_tmp.index:
		if  wl<=right:
			wl_rl=wl
		else:
			wl_rr=wl
			break
	los_rl=series[wl_rl]
	los_rr=series[wl_rr]

	los_l=(left-wl_ll)/(wl_lr-wl_ll)*(los_lr-los_ll)+los_ll
	los_r=(right-wl_rl)/(wl_rr-wl_rl)*(los_rr-los_rl)+los_rl

	# return los_l,los_r
	series=series[series.index<=right]
	series=series[series.index>=left]
	series.loc[left]=round(los_l,3)
	series.loc[right]=round(los_r,3)
	# print(series)
	# tmp=pd.DataFrame([[left,round(los_l,3)],[right,round(los_r,3)]],columns=df.columns)
	# df=df.append(tmp)
	# df=df.set_index('WL')
	# df=df.sort_index()
	# print df
	return series

def APR(file):                #返回分析后的测试数据，包括各通道中心波长及对应插损，通带范围，与其它各通道的串扰，及总中心波长，间隔
	columns=[i+' ITU CT (dB)' for i in list(file.columns)]
	apr=pd.DataFrame(np.zeros((len(columns),len(columns)+6)),
							index=file.columns,
							columns=['PEAK_WL (nm)','IL (dB)','ITU IL (dB)',\
							'1dB BW (nm)','ITU left (nm)','ITU right (nm)']\
							+columns)

	try:
		gap=(find_peak_array(file['CH 1'])[0]-find_peak_array(file['CH 70'])[0])/69
	except:
		gap=(find_peak_array(file['CH 1'])[0]-find_peak_array(file['CH 12'])[0])/11
	result=[]
	for i in file.columns:             #计算通带串扰
		# print type(i) #str
		print(i)      # str CH 1
		ch = int(re.match('CH (\d{1,2})',i).group(1))
		peak_wl,peak_loss=find_peak_array(file[i]) 
		left,right=round(peak_wl-gap/8,3),round(peak_wl+gap/8,3)
		print(left,peak_wl,right)
		# print 
		ITU_upper_min=min(region_interpolation(file[i],left,right))
		tmp=[]
		for j in file.columns:
			if i==j:
				tmp.append('----')
			else:
				ITU_lower_max=max(region_interpolation(file[j],left,right))
				tmp.append(ITU_lower_max-ITU_upper_min)
		# ITU_CT=np.array(tmp)-ITU_upper_min
		bw=bandwidth(file[i],peak_wl,peak_loss,depth=1)
		result.append([peak_wl,peak_loss,ITU_upper_min,bw,left,right]+tmp)
	result=np.array(result)
	result.shape=(len(columns),-1)
	for index in range(len(apr.columns)):
		apr.iloc[:,index]=result[:,index]

	cwl=(float(apr['PEAK_WL (nm)'][0])+float(apr['PEAK_WL (nm)'][-1]))/2

	apr.loc['']=[np.NAN]*len(apr.columns)
	apr.loc['Center WL(nm)']=[cwl]+[np.NAN]*(len(apr.columns)-1)
	apr.loc['Channel gap(nm)']=[round(gap,3)]+[np.NAN]*(len(apr.columns)-1)
	apr.loc['Channel uniformity(dB)']=[apr['IL (dB)'].astype('float').max()-apr['IL (dB)'].astype('float').min()]+[np.NAN]*(len(apr.columns)-1)
	return apr
	# print result
	# print tmp


def get_pure_list(path):                 #返回芯片光谱数据的路径，不包括宽谱光源数据的路径（测试代码时宽谱光源数据名为'ase.csv'，可根据需要修改）
	path_list=os.listdir(path)
	# print path_list
	path_list=filter(lambda x:os.path.isfile(os.path.join(path,x)) and x.split('.')[1].upper()=='CSV' and x!='ase.csv',path_list)
	path_list=sorted(path_list,key=lambda x:int(re.match(r'.*?(\d{1,2}).CSV',x.upper()).group(1)))
	path_list=[os.path.join(path,x) for x in path_list]
	return path_list

def my_mkdir(path):      #创建文件夹
	i=1
	while i!=len(path.split('\\'))+1 and os.path.isdir('//'.join(path.split('\\')[:i])):
		i+=1
	while i!=len(path.split('\\'))+1:
		# print '//'.join(path.split('\\')[:i])
		os.mkdir('//'.join(path.split('\\')[:i]))
		i+=1

if __name__=='__main__':

	# ase=u'G:\\AWG0404\\补测5.9\\NEW.CSV'
	ase=r'D:\a\AWG0404\code\ase.csv'        #纯ASE光谱
	ase=r'F:\小片\ASE.CSV'        #纯ASE光谱
	# path_3ch_origin=r'G:\AWG0404\AWG TMP\237-5'
	path=input('path of file : ')            #芯片光谱存放的文件夹，如D:\a\AWG0404\AWG测试\222-1\ch-1

	path_list=get_pure_list(path)       #单个通道下的所有输出

	try:
		# print path_3ch_origin.split('\\')[-1]
		# print(r'D:\a\AWG0404\data 10.22',path.split('\\')[-2]+'\\'+path.split('\\')[-1])
		new_dir=os.path.join(r'D:\a\AWG0404\data 11.26',path.split('\\')[-2]+'\\'+path.split('\\')[-1])
		my_mkdir(new_dir)
	except Exception as e:
		print(e)
		pass

	acpl=ACPL(path_list,ase)

	path_acpl=os.path.join(new_dir,path.split('\\')[-2]+' '+path.split('\\')[-1]+' ACPL.csv')
	acpl.to_csv(path_acpl)
	

	path_apr=os.path.join(new_dir,path.split('\\')[-2]+' '+path.split('\\')[-1]+' APR.csv')
	APR(acpl).to_csv(path_apr)

	print('\nThe pure spectrum file is in :','\n',path_acpl,'\n')
	print('\nThe analyzed data is saved in :','\n',path_apr,'\n')


