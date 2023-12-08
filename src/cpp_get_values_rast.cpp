#include <Rcpp.h>
#include <math.h>       /* isnan */
#include "windows.h" //Sleep

//' cpp_get_values_rast
//' @name cpp_get_values_rast
//' @param RASTER A RASTER that the value will be extracted from.
//' @param ENVEXT A vector of extent and resolution for environmental matrices: in order: xmin, xmax, ymin, ymax, xres, yres
//' @param XLOCS A vector of x locations.
//' @param YLOCS A vector of y locations.
//' @return A vector of values equal to the length of XLOCS extracted from RASTER.
//' @details A simple C++ function to extract a value from a RASTER
//'   using two vectors describing location. Note that the counting for x and y
//'   begins at zero.
// [[Rcpp::export]]
std::vector<double> cpp_get_values_rast( Rcpp::NumericMatrix RASTER,
                                    std::vector<double> ENVEXT,
                                    std::vector<double> XLOCS,
                                    std::vector<double> YLOCS) {

   // input values
   int nLocs = XLOCS.size();
   int ncols = RASTER.ncol();
   int nrows = RASTER.nrow();
   double xmin = ENVEXT[0];
   double ymin = ENVEXT[2];
   double xres = ENVEXT[4];
   double yres = ENVEXT[5];

   // output vector declaration
   std::vector<double> OUTPUT_VALUES(nLocs);

   for(int loc = 0; loc < nLocs; loc++){

     // transform x,y location to row,col in raster matrix
     int xinx = floor((XLOCS[loc] - xmin) / xres);
     int yinx = floor((YLOCS[loc] - ymin) / yres);
     // std::cout << "  raster location: " << xinx << ", " << yinx << std::endl;

     // CHECK if animal leaves environmental data area
     if( ( xinx > ncols ) |
         ( yinx > nrows ) |
         ( xinx < 0) |
         ( yinx < 0) ){
       // give a value of -99.9 for exceeding environment so the chances of moving
       // are super near zero
       OUTPUT_VALUES[loc] = -99.9;
       Rcpp::Rcerr << "Animal considering exceeding environmental limits, -99.9 value returned instead\n";

     // CHECK for raster NaN
     } else if ( isnan(RASTER(xinx, yinx)) ) {
       OUTPUT_VALUES[loc] = -99.9;
       Rcpp::Rcerr << "Environmental raster returning NaN, -99.9 value returned instead\n";

     // OTHERWISE
     } else {
       // Sleep(100);
       // std::cout << "  raster value: " << RASTER(xinx, yinx) << std::endl;
       OUTPUT_VALUES[loc] = RASTER(xinx, yinx);
     }
   }
   // Sleep(1000);
   return(OUTPUT_VALUES);
 }
