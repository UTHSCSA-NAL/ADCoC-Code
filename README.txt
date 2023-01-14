## ADCoC: Adaptive Distribution Modeling Based Collaborative Clustering for Disentangling Disease Heterogeneity from Neuroimaging Data

## Code for ADCoC

## To reproduce clustering results in the paper, simply run ADCoC.m. 
   The results (adjusted Rand indices) will be stored in ADCoC_Results_ARI.txt.

## To set up the experiments from scartch, run the following scripts:
	gendata.m: Generate the synthetic data, which will be included in sbjFea.mat;
	ADCoC_Coe.m: Form the distribution of the coefficiens, which will be included in SMtr_SC_X_FC_X.mat;
	ADCoC.m: Perform clustering.

## Dependency:
	Y. Li and A. Ngom, "The non-negative matrix factorization toolbox for biological data mining," Source Code for Biology and Medicine, vol. 8, pp. 1-15, 2013.

## If you find the paper and code useful, please consider to cite:
	Hangfan Liu, Michel J. Grothe, Tanweer Rashid, Miguel A. Labrador-Espinosa, Jon B. Toledo, and Mohamad Habes, 
	“ADCoC: Adaptive Distribution Modeling Based Collaborative Clustering for Dis-entangling Disease Heterogeneity from Neuroimaging Data”, 
	IEEE Transactions on Emerging Topics in Computational Intelligence (TETCI), DOI: 10.1109/TETCI.2021.3136587.

