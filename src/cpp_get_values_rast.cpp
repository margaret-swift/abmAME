#include <Rcpp.h>

//' cpp_get_values_rast
//' @name cpp_get_values_rast
//' @param RASTER A RASTER that the value will be extracted from.
//' @param XLOCS A vector of x locations.
//' @param YLOCS A vector of y locations.
//' @return A vector of values equal to the length of XLOCS extracted from
//'   MATRIX.
//' @details A simple C++ function to extract a value from a RASTER
//'   using two vectors describing location. Note that the counting for x and y
//'   begins at zero.

// [[Rcpp::export]]
std::vector<double> cpp_get_values_rast(Rcpp::NumericMatrix RASTER,
                                   std::vector<double> XLOCS,
                                   std::vector<double> YLOCS) {

  int nLocs = XLOCS.size();
  //// ENVIRONMENTAL OBJECTS //////
  //// TODO: UPDATE THIS TO TAKE LATLON instead of RASTER stuff
  int mcols = RASTER.ncol();
  int mrows = RASTER.nrow();
  //////////
  int yOptIndex;
  int xOptIndex;

  //
  std::vector<double> OUTPUT_VALUES(nLocs);

  for(int loc = 0; loc < nLocs; loc++){

    // rounding the locations to correspond to matrix location
    xOptIndex = std::floor(XLOCS[loc]);
    yOptIndex = std::floor(YLOCS[loc]);

    // if animal leaves environmental data area
    if( (xOptIndex > mcols) |
        (yOptIndex > mrows) |
        (xOptIndex < 0) |
        (yOptIndex < 0) ){
      // give a value of -99.9 for exceeding environment so the chances of moving
      // are super near zero
      OUTPUT_VALUES[loc] = -99.9;
      Rcpp::Rcerr << "Animal considering exceeding environmental limits, -99.9 value returned instead\n";
    } else {
      // TODO: make this extract the raster value from the point given
      OUTPUT_VALUES[loc] = RASTER(xOptIndex, yOptIndex);
    }

  }
  return(OUTPUT_VALUES);
}
