#include <Rcpp.h>
#include <math.h>       /* isnan */

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
std::vector<double> cpp_get_values_rast(Rcpp::NumericMatrix RASTER,
                                    std::vector<double> ENVEXT,
                                    std::vector<double> XLOCS,
                                    std::vector<double> YLOCS) {

   // input values
   int nLocs = XLOCS.size();
   double xmin = ENVEXT[0];
   double xmax = ENVEXT[1];
   double ymin = ENVEXT[2];
   double ymax = ENVEXT[3];
   double xres = ENVEXT[4];
   double yres = ENVEXT[5];

   // output vector declaration
   std::vector<double> OUTPUT_VALUES(nLocs);

   // std::cout << "getting raster info for:\n";

   for(int loc = 0; loc < nLocs; loc++){
     // std::cout << "  " << XLOCS[loc] << ", " << YLOCS[loc] << std::endl;

     // if animal leaves environmental data area
     if( ( XLOCS[loc] > xmax ) |
         ( YLOCS[loc] > ymax ) |
         ( XLOCS[loc] < xmin ) |
         ( YLOCS[loc] < ymin ) ){
       // give a value of -99.9 for exceeding environment so the chances of moving
       // are super near zero
       OUTPUT_VALUES[loc] = -99.9;
       Rcpp::Rcerr << "Animal considering exceeding environmental limits, -99.9 value returned instead\n";
     } else {
       // transform xloc and yloc into matrix rows and columns using min val and res
       int xloc = std::floor((XLOCS[loc] - xmin) / xres);
       int yloc = std::floor((YLOCS[loc] - ymin) / yres);

       // grab values from the raster!
       double val = RASTER(xloc, yloc);

       // special check for NaN
       if ( isnan(RASTER(xloc, yloc)) ) {
         val = -99.9;
         Rcpp::Rcerr << "Animal raster returning NaN, -99.9 value returned instead\n";
       }
       OUTPUT_VALUES[loc] = val;
       // std::cout << "  raster val: " << val << std::endl;
     }

   }
   return(OUTPUT_VALUES);
 }
