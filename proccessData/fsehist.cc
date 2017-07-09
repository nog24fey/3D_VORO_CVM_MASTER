#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

int main(int argc, char **argv){

   const string imput1 = "/Users/noguchihironobu/Dropbox/todo0623/FSEngar";
   const string imput2 = "er";
   const string imput3 = "sm1is20000fs1000.dat";

   //vector<string> erate = {"1", "10", "100", "1000"};
   vector<string> erate = {"10", "100"};
   vector<string> area = {"4000", "4100", "4200", "4300", "4400", "4500", "4600", "4700", "4800", "4900", "5000", "5100", "5200", "5270", "5280", "5285", "5290", "5295", "5300", "5305", "5310", "5315", "5320", "5330", "5340", "5350", "5360", "5370", "5380", "5381", "5382", "5383", "5384", "5385", "5386", "5387", "5388", "5389", "5390", "5391", "5392", "5393", "5394", "5395", "5396", "5397", "5398", "5399", "5400", "5410", "5420", "5430", "5450", "5460", "5470", "5480", "5490", "5500", "5600", "5700", "5800", "5900", "6000", "6100", "6200", "6300"};
   for (auto er : erate) {
      for (auto ar : area) {
         double av = 0.0;
         double count = 0.0;
         const string filename = imput1+ar+imput2+er+imput3;
         ifstream ifs;
         ifs.open(filename.c_str());
         //cout<<filename.c_str()<<endl;

         int a1,b1;
         ifs>>a1;
         ifs>>b1;
         int ct = 0;
         while (!ifs.eof()) {
            int a2 = a1;
            int b2 = b1;
            vector<double> arate, e, gbg;
            while ( (a1==a2) && (b1==b2) ) {
               double at,et,gt;
               ifs>>at;
               ifs>>et;
               ifs>>gt;
               //cout<< a1 << " " << b1 << " " << at << " " << et << " "<< gt<<endl;
               ++ct;
               arate.push_back(at);
               e.push_back(et);
               gbg.push_back(gt);
               if ( ifs.eof() ) break;
               ifs>>a2;
               ifs>>b2;
            }
            double e0 = *(e.end()-1);
            if ( e0 < 0.0 ) e0 = 0.0;
            av += e0;
            count += 1.0;
            a1 = a2; b1 = b2;

         }
         ifs.close();
         av /= count;
         cout<< er << " " << ar << " " << av << count << endl;
      }
      cout<<endl;
   }
      
   return 0;	
}
