#include <Rcpp.h>
#include "cpp_get_values_rast.h"
#include "windows.h"
#include <vector>
#include <numeric>

using namespace std;

//' Check if a line segment (step option) crosses a barrier
//' @name check_intersection
//' @param origin start location
//' @param target target location
//' @param BARRIERS  barrier segment endpoint data
//' @param ENVEXT environmental extent for raster data
//' @param LOOKUP lookup table for barrier crossings
//' @return A Boolean of whether any barrier is crossed
//' Useful resource: https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/



// CUSTOM STRUCTURES ===========================================================

// declare a Point structure for p1 and q1
struct Point
{
  double x;
  double y;
};

// END STRUCTURES ==============================================================


// HELPER FUNCTIONS ============================================================
bool onSegment(Point p, Point q, Point r) {
  // Given three collinear points p, q, r, the function checks if
  // point q lies on line segment 'pr'
  if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) &&
     q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))
    return true;
  return false;
}
int orientation(Point p, Point q, Point r) {
  // To find orientation of ordered triplet (p, q, r).
  // The function returns following values
  // 0 --> p, q and r are collinear
  // 1 --> Clockwise
  // 2 --> Counterclockwise
  // See https://www.geeksforgeeks.org/orientation-3-ordered-points/
  // for details of below formula.
  // MAKE SURE THIS IS A DOUBLE AND NOT AN INTEGER!
  double val = ( (q.y - p.y) * (r.x - q.x) ) - ( (q.x - p.x) * (r.y - q.y) );
  if (val == 0) return 0;  // collinear
  return (val > 0)? 1: 2; // clockwise or counterclockwise
}
bool doIntersect(Point p1, Point q1, Point p2, Point q2) {
  // The main function that returns true if line segment 'p1q1'
  // and 'p2q2' intersect.

  // Find the four orientations needed for general and
  // special cases
  int o1 = orientation(p1, q1, p2);
  int o2 = orientation(p1, q1, q2);
  int o3 = orientation(p2, q2, p1);
  int o4 = orientation(p2, q2, q1);

  // General case
  if (o1 != o2 && o3 != o4) return 1;

  // Special Cases
  // p1, q1 and p2 are collinear and p2 lies on segment p1q1
  if (o1 == 0 && onSegment(p1, p2, q1)) return 1;
  // p1, q1 and q2 are collinear and q2 lies on segment p1q1
  if (o2 == 0 && onSegment(p1, q2, q1)) return 1;
  // p2, q2 and p1 are collinear and p1 lies on segment p2q2
  if (o3 == 0 && onSegment(p2, p1, q2)) return 1;
  // p2, q2 and q1 are collinear and q1 lies on segment p2q2
  if (o4 == 0 && onSegment(p2, q1, q2)) return 1;

  return 0; // Doesn't fall in any of the above cases
}
// bool findLookups(double inx,
bool findLookups (int inx,
                  Rcpp::NumericMatrix INXMAT) {
  double n = 10.0;
  double m = 10.0;
  double p = static_cast<double>(INXMAT.ncol());

  // grab leftmost corners of search box
  double start_1 = static_cast<int>( inx - n - (m * p) );
  double start_2 = static_cast<int>( inx - n + (m * p) );

  // if either is negative, set to zero.
  if (start_1 < 0) {start_1 = 0; }
  if (start_2 < 0) {start_2 = 0; }
  std::cout << "left corners: " << start_1 << ", " << start_2 << std::endl;

  // generate sequence of left edge of search box
  // std::vector<double> edge = std::iota( std::begin( start_1 ),
  //                                       std::end( start_2 ),
  //                                       p );

  // for each point on left edge, generate sequence along search box row
  // std::iota( std::begin( xmin ), std::end( xmax ), 0 );
  return true;
}
// END HELPER FUNCTIONS ========================================================



// EXPORTED FUNCTION ===========================================================

// [[Rcpp::export]]
double cpp_check_intersection(std::vector<double> origin,
                              std::vector<double> target,
                              int inx,
                              Rcpp::NumericMatrix BARRIERS,
                              Rcpp::NumericMatrix LOOKUP,
                              Rcpp::NumericMatrix INXMAT
                        ) {

  // declare output and within-function data types
  struct Point p, q, barrier_p, barrier_q;
  p = {origin[0], origin[1]};
  q = {target[0], target[1]};

  // // get raster value for fence-crossing lookup
  // std::cout << "POINT INDEX: " << inx << std::endl;
  // bool answer = findLookups(inx, INXMAT);
  // std::cout << answer << std::endl;

  // calculate a vector of adjacent indices in lookup table
  // std::vector<double> xvals = origin[0];
  // std::vector<double> yvals = origin[0];
  // std::vector<int> v(100) ; // vector with 100 ints.
  // std::iota (std::begin(v), std::end(v), 0); // Fill with 0, 1, ..., 99.

  // use those indices to filter the barrier
  // std::vector<double> ids = BARRIERS(,5);
  // std::cout << "IDS example: " << ids[1] << std::endl;

  // std::vector<int> rows_include = {1,2,3,4};// find which ids are in our vector created above
  // std::vector<int> rows_include = {1,2,3,4};
  // Rcpp::NumericMatrix BARRIERS_FILT = BARRIERS(rows_include,);

  // std::cout << "SUB MATRIX: " << BARRIERS_FILT << std::endl;



  // Loop over barrier segments and see if there's a crossing
  int L = BARRIERS.nrow();
  for (int i = 0; i < L; i++) {

    // get info for this barrier segment
    barrier_p = { BARRIERS(i,0), BARRIERS(i,1) }; //{ barrier_x1[i], barrier_y1[i] };
    barrier_q = { BARRIERS(i,2), BARRIERS(i,3) }; //{ barrier_x2[i], barrier_y2[i] };

    // Animal should step through the barrier with P=p_cross. Draw a random number
    // from the uniform distribution and compare to p_cross for this barrier segment
    bool doesIntersect = doIntersect(p, q, barrier_p, barrier_q);
    double p_cross = BARRIERS(i,4);
    double d_cross = Rcpp::runif(1, 0, 1)[0];

    // if d_cross is ABOVE the p_cross threshold, THROW AWAY current draw.
    if ( doesIntersect && d_cross > p_cross ) { return true; }
  }

  // If there are no crossings, then individual is not blocked
  return false;
}
// END EXPORTED FUNCTION =======================================================

