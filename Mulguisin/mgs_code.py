from numba import njit
from numba.typed import List

@njit
@njit
def mgs3d(Rcut,x1,y1,z1,isort):
    
	t_start = time()
    
	Nmgs=0
	imgs = List()
	cng = List()
	clg = List()
	clm = List()
	#cll = List()
	for i in range(len(x1)):
		ic = isort[i]
		Icc=-1
		Iwg=-1
		# Dmin larger than Rcut**2
		Dmin=Rcut**2 + 1

		for ig in range(0,Nmgs):
			Dist=Rcut**2 + 1
			for icl in range(0,cng[ig]):
				itmp=clg[ig][icl]
				d2 = (x1[ic]-x1[itmp])**2+(y1[ic]-y1[itmp])**2+(z1[ic]-z1[itmp])**2
				if (d2<Dist): 
					idmin=itmp
					Dist=d2
			if (Dist < Dmin):
				Iwg=ig
				Icc=idmin
				Dmin=Dist
		if (Dmin < Rcut*Rcut):
			imgs.append(Iwg)
			clg[Iwg].append(ic)
			clm[Iwg].append(Icc)
			cng[Iwg]+=1
		else:
			imgs.append(Nmgs)	
			clg.append([ic])
			clm.append([ic])
			cng.append(1)
			Nmgs+=1
    
	t_end = time()
	print('Calculation is done. Time = ',t_end - t_start)
    
	return Nmgs,imgs,clg,clm,cng

@njit
def mgs3d_period(Rcut,x1,y1,z1,isort,Lbox):
    
#	t_start = time()
    
	Nmgs=0
	imgs = List()
	cng = List()
	clg = List()
	clm = List()
	#cll = List()
	for i in range(len(x1)):
		ic = isort[i]
		Icc=-1
		Iwg=-1
		# Dmin larger than Rcut**2
		Dmin=Rcut**2 + 1

		for ig in range(0,Nmgs):
			Dist=Rcut**2 + 1
			for icl in range(0,cng[ig]):
				itmp=clg[ig][icl]
                		#We are calculating the distance, so only the absolute value is necessary.
				del_x = np.abs(x1[ic]-x1[itmp])
				del_y = np.abs(y1[ic]-y1[itmp])
				del_z = np.abs(z1[ic]-z1[itmp])
                
                		#Periodic boundary condition. We assume that the cluster size is smaller than 0.5*Lbox
				if del_x > Lbox*0.5:
				    del_x = Lbox - del_x
				if del_y > Lbox*0.5:
				    del_y = Lbox - del_y
				if del_z > Lbox*0.5:
				    del_z = Lbox - del_z
    
                
				d2 = (del_x)**2+(del_y)**2+(del_z)**2                
        
				if (d2<Dist): 
					idmin=itmp
					Dist=d2
			if (Dist < Dmin):
				Iwg=ig
				Icc=idmin
				Dmin=Dist
		if (Dmin < Rcut*Rcut):
			imgs.append(Iwg)
			clg[Iwg].append(ic)
			clm[Iwg].append(Icc)
			cng[Iwg]+=1
		else:
			imgs.append(Nmgs)	
			clg.append([ic])
			clm.append([ic])
			cng.append(1)
			Nmgs+=1
	return Nmgs,imgs,clg,clm,cng
