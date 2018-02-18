import numpy as np
import scipy as sp
import matplotlib as plt
#construction of the network and its elements
class Network:
	#do we need to mark reactants and products??????????
	class Cell:#initial state
		def __init__(self):
			self.const = np.array([0,0])#random kinetic constant,eqauls 0 initially:theta,gamma
			self.proteinlist = []
			for i in range(2):
				proteinlist.append(protein())
			self.proteincomplex = np.array([])
			self.score = 0;
			self.amount = 4#intial = 4 concentration = 1/4 each 




		class protein:# gene and protein are 1-1 initially
			class gene:
				pass
			
			def __init__(self):
				self.rates = []##tau,delta
				rates.append((-1)*np.random.uniform(0,1))
				rates.append(np.random.uniform(0,1))
				self.gene = []
				gene.append(gene())

		def mutateinCell(self,n):
		#5 cases:
			situation = randi(1,5)
			if (situation == 1):
				changekinetic()
			else if (situation == 2):
				changeDegradeorProduction()
			else if (situation == 3):
				addGene()
			else if (situation ==4):##3 reactions happen silutaniously or choose 1????
				c = randint(0,1)
				proteinlist.append(protein())

				if (c==0):#not reverse
					index1 = np.random.randint(0,len(proteinlist)-1)
					index2 = np.random.randint(0,len(proteinlist)-1)
					p1 = proteinlist[index1]
					p2 = proteinlist[index2]#find another gene promoter,which is a protein
					del proteinlist[index1]
					del proteinlist[index2]
					proteincomplex.append([p1,p2])
					p2.rates[0] = np.random.rand(0,1)
				else:#reverse
					index = proteincomplex[np.random.randint(0,len(proteincomplex)-1)]
					proteinlist.append(proteincomplex[0])
					proteinlist.append(proteincomplex[1])
					np.delete(proteincomplex,index)#?????????do we need to delete???????????????????????/
			else:#situation5
				singleprotein = np.random.randint(0,1)#choose whether it's one or two case
				if (singleprotein == 1):
					single=np.random.randint(0,1)#divides into 4 situations
					if (single==0):#single protein;modify tau,adding protein and gene for each cell
						addProtein()

					else :#protein complex;choose 1 complex,decompose,add survival one into proteinlist
						modify=np.random.randint(0,1)
						idx=np.random.randint(0,proteincomplex.size-1)
						proteinlist.append(proteincomplex[idx][modify])
						np.delete(proteincomplex,idx)

				else:#2 protein case---what if 2 proteincomplexes are chosen?????
					size = len(proteinlist)+proteinlist.size
					index1 = np.random.randint(0,size)#first element choosen
					index2 = np.random.randint(0,size)#second element choosen
					while (index1>len(proteinlist)&&index2>len(proteinlist)):
						index1 = np.random.randint(0,size)
						index2 = np.random.randint(0,size)				
					if (index1<proteincomplex.size()&&index1>len(proteinlist))||(index2<proteincomplex.size()&&index2>len(proteinlist)):#choose 2 proteins
						#only 1 is protein
						if index1<proteincomplex.size()&&index1>len(proteinlist):
							survivor = np.random.sample(proteincomplex[index1],1)
							survivor = np.random.sample(survivor.append(proteinlist[index2]),1)
							const[0] = np.random.randint(0,1)
						else:
							survivor = np.random.sample(proteincomplex[index2],1)
							survivor = np.random.sample(survivor.append(proteinlist[index1]),1)
							const[0] = np.random.randint(0,1)						
					else :
						#choose 2 protein complex or 2 proteins
						c = np.random.randint(0,1)#dimerization or degradation
						if c==0:
							#dimerization
							const[1]=np.random.rand(0,1)
							proteincomplex.append([proteinlist[index1],proteinlist[index2]])
						else:
							survivor = (c==1)?proteinlist[index2]:proteinlist[index1]
							const[1] = np.random.rand(0,1)

	
	def __init__(self,n):
		self.cells = []
		self.size = n
		for i in range(n):
			elem = Cell()
			self.cells.append(elem)


	def changekinetic(self):
		for i in cells:
			one = randint(0,3)
			if one==1 or one==0:
				for j in i.proteinlist:
					if one==1:
						j.rates[1] *= np.random.uniform(0,2)
					else:
						j.rates[0] *= np.random.uniform(0,2)
			else:
				if one==2:
					if i.const[0]==0:
						i.const[0] = np.random.rand(0,1)
					else:
						i.const[0] *= np.random.uniform(0,2)
				else:
					if i.const[1]==0:
						i.const[1] = np.random.rand(0,1)
					else:
						i.const[1] *= np.random.uniform(0,2)					



	def changeDegradeorProduction(self):
		for i in cells:
			for j in i.proteinlist:
				j.rates[0] *= np.random.uniform(0,2)
				j.rates[1] *= np.random.uniform(0,2)

	def addGene(self):##does gene and protein needs to occur pairly?suppose add.
		for i in cells:
			i.proteinlist.append(protein())
			amount+=2

	def addProtein(self):
		for i in cells:
			i.proteinlist.append(protein())
			amount+=2

#mutation part:

	def mutate(self):
		copy = Network(self.size)
		for i in copy.cells:
			i.mutateinCell()



	def calculateScore(self,cell):
		rv = 0




	
