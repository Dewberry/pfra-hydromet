import unittest
import sys
import pandas as pd
sys.path.append('../core')
from core import hydromet 
print("heloo")

class TestHydroMet(unittest.TestCase):
  
    @classmethod
    def setUpClass(cls):
        # print('setupClass')   
        pass

    @classmethod
    def tearDownClass(cls): 
        # print('teardownClass')
        pass

    def setUp(self):
        '''Sample data for test case'''
        self.Tr = [2,5,10,25,50,100,200,500,1000]
        self.data = {'Lower (90%)':    [2.43, 3.03, 3.56, 4.30, 4.86, 5.39, 5.91, 6.71, 7.32],
                     'Expected Value': [2.83, 3.55, 4.20, 5.20, 6.07, 7.03, 8.09, 9.62, 10.89],
                     'Upper (90%)':    [3.33, 4.19, 4.98, 6.50, 7.66, 9.07, 10.71, 13.07, 14.85]}
        self.df = pd.DataFrame(self.data, index=self.Tr)
        self.df.index.name='Tr'
        self.seed = 100

    def tearDown(self):
        # print('tearDown')
        pass
    
    def test_extrap_add_ari(self):
        print('Testing','test_extrap_add_ari')

        df_new = hydromet.extrap_add_ari(self.df,display_print=False)

        for idx in self.Tr :
            for col in self.df.columns:
                self.assertEqual(self.df.loc[idx,col], df_new.loc[idx,col])

        for col in self.data.keys():
            self.assertGreaterEqual(df_new.loc[3000,col], self.df.loc[1000,col])


    def test_add_ari(self):
        print('Testing','add_ari')
        df_new = hydromet.add_ari(self.df)
        check_cols = ['Ann. Exc. Prob.', 'ARI', 'Log10_ARI']
        
        for idx in self.Tr :
            for col in check_cols:
                self.assertGreaterEqual(df_new.loc[idx,col], self.df.loc[idx,col])


      
if __name__ == '__main__':

    unittest.main()