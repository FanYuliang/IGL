import numpy as np
import scipy as sp
import matplotlib as plt
#construction of the network and its elements
class Network:
	def _self_(self,n):
		self.cells = []
		self.size = n
		for i in range(n):
			elem = Cell()
			self.cells.append(elem)


	def changekinetic(self):
		for i in cells:
			one = randint(0,3)
			i.const[one] *= np.random.uniform(0,2)



	def changeDegradeorProduction(self):
		for i in cells:
			one = randint(0,3)
			i.p2.degrade *= np.random.uniform(0,2)
			i.p1.degrade *= np.random.uniform(0,2)

	def addGene(self):##does gene and protein needs to occur pairly?suppose add.
		for i in cells:
			i.genelist.append(gene())
			i.proteinlist.append(protein())



	class Cell:#initial state
		def _self_(self):
			self.const = np.random.rand(1,4)#random kinetic constant,tau,delta,theta,gamma
			self.proteinlist = []
			self.genelist = []
			for i in range(2):
				proteinlist.append(protein())
				genelist.append(gene)
			self.complex = []
			self.amount = 4#intial = 4 concentration = 1/4 each 

		class gene:#nothing but to denote the gene element
			pass


		class protein:
			def _self_(self,complex):#denotes whether it's a protein complex
				self.product = np.random.rand()##?????do we need product rate?
				self.degrade = np.random.rand()

#mutation part:
	def mutate(self,n):
		#5 cases:
		copy=Network(n.size)
		situation = randi(1,5)
		if (situation == 1):
			changekinetic()
		else if (situation == 2):
			changeDegradeorProduction()
		else if (situation == 3):
			addGene()
		else if (situation ==4):##3 reactions happen silutaniously or choose 1????
			c = randint(0,1)
			newgene = gene()
			newprotein = protein()
			proteinlist.append(newprotein)
			genelist.append(newgene)

			if (c==0):
				p = proteinlist[randint(0,proteinlist.size()-1)]
				g = genelist[randint(0,proteinlist.size()-1)]

			else:
				complex1 = complex[randint(0,complex.size()-1)]
				proteinlist.append(complex1[0])
				proteinlist.append(complex1[1])
				complex.remove(complex1)


	def calculateScore(self,prot,):



	
