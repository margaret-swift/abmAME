#include <Rcpp.h>
#include <cmath>
#include "windows.h" // for sleep, i.e. 2 secs = Sleep(2);
#include "cpp_vonmises.h"
#include "cpp_sample_options.h"
#include "cpp_cycle_draw.h"
#include "cpp_get_values_rast.h"
#include "cpp_maxmin.h"
#include "cpp_check_intersection.h"

 //' @name cpp_abm_simulate
 //' @title cpp_abm_simulate
 //' @description The C++ function that runs the agent based animal movement
 //'   model.
 //' @param startx The x coordinate start location
 //' @param starty The y coordinate start location
 //' @param timesteps The number of timesteps to be simulated
 //' @param ndes The number of possible destination considered for foraging.
 //' @param nopt The number of options to be considered at each timestep
 //' @param shelter_locs_x A vector of shelter locations x coordinates.
 //' @param shelter_locs_y A vector of shelter locations y coordinates.
 //' @param sSiteSize A double dictating at what range the animal's movements
 //'   will dramatically drop simulating (near-)stationary behaviour.
 //' @param avoidPoints_x A vector of avoidance point x coordinates.
 //' @param avoidPoints_y A vector of avoidance point y coordinates.
 //' @param k_desRange Shape parameter describing the gamma distribution
 //'   destinations can be draw from.
 //' @param s_desRange Scale parameter describing the gamma distribution
 //'   destinations can be draw from.
 //' @param mu_desDir Mean of the turn angle destinations are drawn from.
 //' @param k_desDir Concentration of the turn angles destinations are drawn
 //'   from.
 //' @param destinationTrans 0 - no transformation applied to the distance to
 //'   destination weighting, 1 - distance to destination weighing is
 //'   square-rooted, 2 - distance to destination weighting is squared
 //' @param destinationMod A coefficient to be applied to the distance to
 //'   destination weighting.
 //' @param avoidTrans 0 - no transformation applied to the distance to avoidance
 //'   points weighting, 1 - distance to avoidance points weighing is
 //'   square-rooted, 2 - distance to avoidance points weighting is squared
 //' @param avoidMod A coefficient to be applied to the avoidance points
 //'   weighting.
 //' @param k_step Parameter describing step length
 //' @param s_step Parameter describing step
 //' @param mu_angle Parameter describing angle
 //' @param k_angle Parameter describing angle variation
 //' @param rescale A value that describes the cell size of the environmental
 //'   matrices relative to the units of the step lengths. Must be greater than
 //'   zero. The step lengths returned are on the scale of the matrix and need to
 //'   be back transformed to match the input step length units for comparison.
 //'   Default is 1, step and matrix unit are the same.
 //' @param b0_Options Behave transitional probs for behave 0
 //' @param b1_Options Behave transitional probs for behave 1
 //' @param b2_Options Behave transitional probs for behave 2
 //' @param rest_Cycle_A The amplitude of the resting/active cycle.
 //' @param rest_Cycle_M The cycle offset from 0 (Midline Statistic Of Rhythm) of
 //'   the resting/active cycle.
 //' @param rest_Cycle_PHI The offset of the cycle (\eqn{\phi}; i.e., acrophase)
 //'   of the resting/active cycle.
 //' @param rest_Cycle_TAU Cycle frequency (\eqn{\tau}; i.e., period) of the
 //'   resting/active cycle.
 //' @param addCycles A integer describing how many additional activity cycles to
 //'   include, can be zero.
 //' @param add_Cycle_A The amplitude of the additional active cycles, should be
 //'   a vector of doubles.
 //' @param add_Cycle_M The cycle offset from 0 (Midline Statistic Of Rhythm) of
 //'   the additional active cycles, should be a vector of doubles.
 //' @param add_Cycle_PHI The offset of the cycle (\eqn{\phi}; i.e., acrophase)
 //'   of the additional active cycles, should be a vector of doubles.
 //' @param add_Cycle_TAU Cycle frequency (\eqn{\tau}; i.e., period) of the
 //'   additional active cycles, should be a vector of doubles.
 //'
 //' @param shelterMatrix The RASTER describing shelter site quality.
 //' @param forageMatrix The RASTER describing foraging site quality.
 //' @param moveMatrix The RASTER describing movement ease.
 //'
 //'---------------------------------------------------
 //' MAGGIE'S ADDITIONAL PARAMETERS FOR FENCE PROJECT
 //' @param envExt extent and resolution for environmental matrices;
 //'                      in order: xmin, xmax, ymin, ymax, xres, yres
 //' @param barriers = barriers table with columns: x, y, xend, yend, perm, id
 //'---------------------------------------------------
 //'
 //' @return A list of simulated animal details to be passed and handled by R
 //'   function to make easier to use.
 //' @details Requires cmath and headers for smaller C++ functions for draws.
 //' @seealso [abm_simulate()] is the R function wrap with more in depth
 //'   documentation.

 // [[Rcpp::export]]
 Rcpp::List cpp_abm_simulate(
     double startx,
     double starty,
     int timesteps,
     int ndes,
     int nopt,

     std::vector<double> shelter_locs_x,
     std::vector<double> shelter_locs_y,
     double sSiteSize,
     std::vector<double> avoidPoints_x,
     std::vector<double> avoidPoints_y,

     double k_desRange,
     double s_desRange,
     double mu_desDir,
     double k_desDir,
     int destinationTrans,
     double destinationMod,
     int avoidTrans,
     double avoidMod,

     std::vector<double> k_step,
     std::vector<double> s_step,
     std::vector<double> mu_angle,
     std::vector<double> k_angle,
     double rescale,

     std::vector<double> b0_Options,
     std::vector<double> b1_Options,
     std::vector<double> b2_Options,

     double rest_Cycle_A,
     double rest_Cycle_M,
     double rest_Cycle_PHI,
     double rest_Cycle_TAU,

     int addCycles,
     std::vector<double> add_Cycle_A,
     std::vector<double> add_Cycle_M,
     std::vector<double> add_Cycle_PHI,
     std::vector<double> add_Cycle_TAU,

     Rcpp::NumericMatrix shelterMatrix,
     Rcpp::NumericMatrix forageMatrix,
     Rcpp::NumericMatrix moveMatrix,
     std::vector<double> envExt,
     Rcpp::NumericMatrix barriers
 ){

   // DISTRIBUTION DRAW OBJECTS -------------------------------------------------
   // location for the step and angle for each choice or destination
   double step;
   double vmdraw;
   double angle;
   //----------------------------------------------------------------------------

   // MOVEMENT OPTIONS OBJECTS --------------------------------------------------
   // the options stores for the loop
   std::vector<double> x_Options(nopt);
   std::vector<double> y_Options(nopt);
   std::vector<int> step_Options(nopt);
   std::vector<double> move_Options(nopt);
   // needed for the choosing of option
   int chosen;
   // store the choice at each step
   std::vector<int> chosen_Options(timesteps);
   // the options stores for the output including all options
   std::vector<double> x_OptionsAll(nopt*timesteps);
   std::vector<double> y_OptionsAll(nopt*timesteps);
   std::vector<int> step_OptionsAll(nopt*timesteps);
   std::vector<double> stepAll(nopt*timesteps);
   // the destination stores for the output including all destination, they will but much longer than needed
   std::vector<double> x_DesOptionsAll(nopt*timesteps);
   std::vector<double> y_DesOptionsAll(nopt*timesteps);
   std::vector<int> step_DesOptionsAll(nopt*timesteps);
   std::vector<int> behave_DesOptionsAll(nopt*timesteps);
   //----------------------------------------------------------------------------

   // REALISED/CHOSEN MOVEMENT OBJECTS ------------------------------------------
   std::vector<double> x_Locations(timesteps);
   std::vector<double> y_Locations(timesteps);
   std::vector<double> sl_Locations(timesteps);
   std::vector<double> ta_Locations(timesteps);
   std::vector<int> step_Locations(timesteps);
   // also recording the destination aimed for at each timestep
   std::vector<double> des_x_Locations(timesteps);
   std::vector<double> des_y_Locations(timesteps);
   std::vector<double> des_chosen(timesteps);
   //----------------------------------------------------------------------------

   // BEHAVIOUR RELATED OBJECTS -------------------------------------------------
   // somewhere to store the behaviours at each step
   std::vector<int> behave_Locations(timesteps);
   // initial behaviour set to 0
   behave_Locations[0] = 0;
   // and the behaviour related movement characteristics for initial state
   double behave_k_step = k_step[0];
   double behave_s_step = s_step[0];
   double behave_mu_angle = mu_angle[0];
   double behave_k_angle = k_angle[0];
   std::string behave_cur = ""; // MAGGIE

   // something that store the time adjusted behavioural shifts
   // and initialise them with the provided base values
   std::vector<double> b0_Options_Current = b0_Options;
   std::vector<double> b1_Options_Current = b1_Options;
   std::vector<double> b2_Options_Current = b2_Options;
   // a store for the modifiers that will change throughout the cycles
   double b0_dailyMod;
   double b0_addMod;
   //----------------------------------------------------------------------------

   // DESINATION OBJECTS -------------------------------------------------------
   int chosenDes = 0;

   // desMatrix will update depending on behaviour
   Rcpp::NumericMatrix desMatrix;


   // define the memory for the shelter site selection
   // requires predefined vector of shelter coords
   int shel_ndes = shelter_locs_x.size();
   std::vector<double> des_Options(shel_ndes);
   // set up the space for foraging destinations that can be dynamically selected
   // IE do not have to be predefined like the shelter ones
   std::vector<double> x_forageOptions(ndes);
   std::vector<double> y_forageOptions(ndes);
   std::vector<double> des_forageOptions(ndes);
   //----------------------------------------------------------------------------

   // AVOIDANCE OBJECTS ---------------------------------------------------------
   double navp = avoidPoints_x.size();
   double cumulative_dist;
   std::vector<double> distance_toAvoid(nopt);
   //----------------------------------------------------------------------------

   // DISTANCE CALCULATION OBJECTS ----------------------------------------------
   // objects for guiding animal towards destination
   double currDist;
   // a vector to hold the distance between options and the destination
   double c_dist2;
   double c_dist;
   std::vector<double> distance_toDes(nopt);
   // double distInvert;
   std::vector<double> weights_toDes(nopt);
   // we need storage for the last angle followed as turning angle is relative
   double last_angle = 0;
   // as we don't know which one is chosen until after selection we need to store
   // all until then
   std::vector<double> taOptions(nopt);
   // and we can add in step to make review easier
   std::vector<double> slOptions(nopt);
   // to initialise the animal is attracted to the first shelter site
   double des_x = shelter_locs_x[0];
   double des_y = shelter_locs_y[0];
   //----------------------------------------------------------------------------

   // INITIAL LOCATION SETTING --------------------------------------------------
   x_Locations[0] = startx;
   y_Locations[0] = starty;
   // intial destination
   des_x_Locations[0] = des_x;
   des_y_Locations[0] = des_y;
   //----------------------------------------------------------------------------
   /*initial options are populated with start location as animal doesn't make an
    initial choice, that way the option indexing matches the timestep value
    within the loop */
   for(int u = 0; u < nopt; u++){
     x_OptionsAll[u] = startx;
     y_OptionsAll[u] = starty;
     step_OptionsAll[u] = 0;
     stepAll[u] = 0;
   }


   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // LOOP OVER EACH TIMESTEP
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   // LOGGING!
   int nstep = 50;
   int stepcount = 0;
   int slice = std::floor(timesteps / nstep);
   bool dolog = timesteps > 1000;
   if (dolog) {
     std::cout << "[" << std::string(nstep, '.') << "]" << std::endl;
     std::cout << " ";
   }

   for(int i = 1, a = nopt, desi = 0; i < timesteps; i++){
     // std::cout << "LOOP " << i << std::endl;

     // progress marker
     if ( (dolog) && (i > 0) && (i % slice == 0)) {
       std::cout << '|';
       stepcount++;
     }

     // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     // CYCLES
     // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     /* working under the assumption that i == minute, but the cycle is defined in
      hours AKA 12 hour cycle offset to be crepusclar, we need to convert i AKA minute to hours */
     b0_dailyMod = cpp_cycle_draw(
       i*1.0 / 60, // make i a double and convert it to hours. i == 1 min so 1/60 i == hour
       rest_Cycle_A,
       rest_Cycle_M,
       rest_Cycle_PHI / rest_Cycle_TAU, // make sure PHI is kept ~ to TAU so no drift
       rest_Cycle_TAU);


     if(addCycles > 0){

       for(int cyc = 0; cyc < addCycles; cyc++){

         b0_addMod = cpp_cycle_draw(
           i*1.0 / 60, // make i a double and convert it to hours.
           // i == 1 min so 1/60 i == hour
           add_Cycle_A[cyc],
                      add_Cycle_M[cyc],
                                 add_Cycle_PHI[cyc] / add_Cycle_TAU[cyc], // make sure PHI is kept ~ to
                                                                   // TAU so no drift
                                                                   add_Cycle_TAU[cyc]);
         // we then update the resting chance modifier with the second cycle output
         b0_dailyMod = b0_dailyMod + b0_addMod;
       }


     } // end of if




     // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     // CHOOSE BEHAVIORAL STATE AND TRANSITION PROB
     // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     /* switch to use a given set of transition probabilities that change
      depending on the previous behavioural state*/
     switch(behave_Locations[i-1]){
     case 0:
       // this will update the behaviour shift prob depending on the time of day
       b0_Options_Current[0] = b0_Options[0] + b0_dailyMod;
       // draw from the updated behaviour probs to get the next behavioural state
       behave_Locations[i] = cpp_sample_options(b0_Options_Current);
       break;
     case 1:
       b1_Options_Current[0] = b1_Options[0] + b0_dailyMod;
       behave_Locations[i] = cpp_sample_options(b1_Options_Current);
       break;
     case 2:
       b2_Options_Current[0] = b2_Options[0] + b0_dailyMod;
       behave_Locations[i] = cpp_sample_options(b2_Options_Current);
       break;
     } // end of switch



     // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     // STEP, ANGLE, AND RASTER MATRIX BASED ON BEHAVIOR
     // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     /* assigning the step and angle parameters
      depending on the behaviour */
     switch(behave_Locations[i]){
     case 0:
       behave_k_step = k_step[0];
       behave_s_step = s_step[0];
       behave_mu_angle = mu_angle[0];
       behave_k_angle = k_angle[0];
       desMatrix = shelterMatrix;
       behave_cur = "sheltering";
       break;
     case 1:
       behave_k_step = k_step[1];
       behave_s_step = s_step[1];
       behave_mu_angle = mu_angle[1];
       behave_k_angle = k_angle[1];
       desMatrix = forageMatrix;
       behave_cur = "foraging";
       break;
     case 2: // MAGGIE note: swapped 1 and 2
       behave_k_step = k_step[2];
       behave_s_step = s_step[2];
       behave_mu_angle = mu_angle[2];
       behave_k_angle = k_angle[2];
       desMatrix = moveMatrix;
       behave_cur = "moving";
       break;
     }


     // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     // CHOOSING APPROPRIATE DESTINATION
     // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     // choosing between a number of possible predefined destinations
     // new destination should be chosen if behaviour changes or after a
     // period of time passes e.g. a day.
     if(!(behave_Locations[i-1] == behave_Locations[i])){

       // weight up options based on quality of location from raster
       // desMatrix updated at each behaviour switch above
       // switch to update possible destination / point of attraction
       switch(behave_Locations[i]){
       case 0:

         // MAGGIE: RASTER CHECK #1
         // std::cout << "resting raster check" << std::endl;
         // Sleep(2);
         des_Options = cpp_get_values_rast(desMatrix, envExt, shelter_locs_x, shelter_locs_y);
         chosenDes = cpp_sample_options(des_Options);

         des_x = shelter_locs_x[chosenDes];
         des_y = shelter_locs_y[chosenDes];

         for(int shelopt = 0; shelopt < shel_ndes; shelopt++, desi++){
           // iterates to next one for following allocation of the destinations chosen from
           x_DesOptionsAll[desi] = shelter_locs_x[shelopt];
           y_DesOptionsAll[desi] = shelter_locs_y[shelopt];
           step_DesOptionsAll[desi] = i;
           behave_DesOptionsAll[desi] = behave_Locations[i];
         }

         break;
       // MAGGIE note: swapped cases 1 and 2 for forage/explore
       case 1:

         for(int dopt = 0; dopt < ndes; dopt++, desi++){

           step = Rcpp::rgamma(1, k_desRange, s_desRange)[0];
           step = step / rescale;
           vmdraw = cpp_vonmises(1, mu_desDir, k_desDir)[0];
           angle = vmdraw * 180/M_PI;
           angle = last_angle + angle;
           x_forageOptions[dopt] = x_Locations[i-1] + cos(angle) * step;
           y_forageOptions[dopt] = y_Locations[i-1] + sin(angle) * step;

           x_DesOptionsAll[desi] = x_forageOptions[dopt];
           y_DesOptionsAll[desi] = y_forageOptions[dopt];
           behave_DesOptionsAll[desi] = behave_Locations[i];
           step_DesOptionsAll[desi] = i;
         }

         // MAGGIE: RASTER CHECK #2
         // std::cout << "forage raster check" << std::endl;
         // Sleep(2);
         des_forageOptions = cpp_get_values_rast(desMatrix, envExt, x_forageOptions, y_forageOptions);

         chosenDes = cpp_sample_options(des_forageOptions);
         des_x = x_forageOptions[chosenDes];
         des_y = y_forageOptions[chosenDes];

         break;

       case 2:
         // we can update the destination here, but explore doesn't have a
         // destination weighting so this has no effect
         des_x = shelter_locs_x[0];
         des_y = shelter_locs_y[0];
         break;
       }
     }

     // store the destination for output
     des_x_Locations[i] = des_x;
     des_y_Locations[i] = des_y;
     des_chosen[i] = chosenDes;


     // current distance from destination, needed to see if movement should
     // reduce because of proximity
     c_dist2 = std::pow((des_x - x_Locations[i-1]), 2) +
       std::pow((des_y - y_Locations[i-1]), 2);
     currDist = std::sqrt(c_dist2);



     //     ***************************************************************
     //   *******************************************************************
     // ***********************************************************************
     //                           STEP OPTIONS DRAW !
     // ***********************************************************************
     //   *******************************************************************
     //     ***************************************************************

     // MAGGIE BOOKMARK : Changed original FOR loop to a WHILE loop so that if
     //   a location is chosen that puts the animal within or on the other side
     //   of a barrier, that option is skipped (with a certain barrier-permeability
     //   probability P) and another is chosen.
     // Step 0: Restrict barrier list to those within Xkm of origin point, to save time(?)
     // TODO ^^ BOOKMARK ^^

     int j = 0;
     int trial_count = 0;
     int trial_kill = nopt * 10000;

     while (j < nopt) {

       // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       // INITIALIZE MOVEMENT PARAMETERS FOR THIS DRAW
       // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       // std::cout << "OPTION J = " << j << " OF " << nopt << std::endl;
       if(j == 0){ // repeat for each start of each step
         /* for each step set the location as the previously chosen location */
         x_Options[0] = x_Locations[i-1];
         y_Options[0] = y_Locations[i-1];
         step_Options[0] = i;
         slOptions[0] = 0;
         taOptions[0] = last_angle;

         // these need assignment regardless
         x_OptionsAll[a] = x_Options[j];
         y_OptionsAll[a] = y_Options[j];
         step_OptionsAll[a] = i;
         j++; a++;// increment before 'continue' or it will run away
         continue;
       }

       // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       // DRAW STEP LENGTH AND TURNING ANGLE
       // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       // an if to make sure the animal doesn't move too far from a shelter site
       // need the initial high timesteps to get there.
       // this could be swapped to maximise the step length if it is far from shelter/centre
       if( (behave_Locations[i] == 0) & (currDist < (sSiteSize/rescale)) ){
         step = Rcpp::rgamma(1, behave_k_step, behave_s_step)[0];
         step = step / rescale / 100;
       } else{
         step = Rcpp::rgamma(1, behave_k_step, behave_s_step)[0];
         step = step / rescale;
       }

       // draw from vonmises for angle
       vmdraw = cpp_vonmises(1, behave_mu_angle, behave_k_angle)[0];
       angle = vmdraw * 180/M_PI;
       angle = last_angle + angle;

       // turning angle and step length options
       taOptions[j] = angle;
       slOptions[j] = step;

       // choose the corresponding x and y coordinates
       x_Options[j] = x_Options[0] + cos(angle) * step;
       y_Options[j] = y_Options[0] + sin(angle) * step;

       // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       // CHECK POTENTIAL BARRIER CROSSING BEHAVIOR
       // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       // step 1: create a line between current point and new point
       double x0 = x_Options[0];
       double y0 = y_Options[0];
       double xj = x_Options[j];
       double yj = y_Options[j];

       // step 2: check if the origin-target line intersects with any barrier nearby
       // https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
       std::vector<double> origin = {x0, y0};
       std::vector<double> target = {xj, yj};

       bool isblocked = cpp_check_intersection( origin, target, // these are individual points
                                                barriers);

       // step 3: if there is an intersection, discard with probability 1-P
       if (isblocked) {
         if (dolog) {
           std::cout << "[" << std::string(nstep, '.') << "]" << std::endl;
           std::cout << " " << std::string(stepcount, '|');
         }

         // If there have been too many trials, keep the most recent draw
         //   and alert the user.
         if (trial_count >= trial_kill) {
           std::cout << "WARNING ele cannot find target that doesn't cross " <<
             "the barrier! Aborting selection at DRAW " << i << ", option " <<
               j << ".\n";
         // otherwise, the intersection blocks the individual; redraw but
         // increase the trial counter
         } else {
           trial_count++;
           continue; // THROW AWAY current draw
         }
       }

       // add in which step the options are for
       step_Options[j] = i;

       // a is keeping tracking of the position in a vector timesteps*nopts
       x_OptionsAll[a] = x_Options[j];
       y_OptionsAll[a] = y_Options[j];
       step_OptionsAll[a] = i;
       stepAll[a] = step;

       // choice vector is needed for the sample function later on
       // choicesVec[j] = j;
       j++;
       a++;
     } // END MOVEMENT LOOP

     // WEIGHTING OPTIONS BY DISTANCE FROM DESTINATION
     // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     // MAGGIE: RASTER CHECK #3
     // std::cout << "fence crossings: " << crosscount << std::endl;
     // std::cout << "checking move raster values at options" << std::endl;
     move_Options = cpp_get_values_rast(moveMatrix, envExt, x_Options, y_Options);

     /* here we need to adjust the movement objects so the animal prefers to head
      * towards the centre point.
      * a2 + b2 = c2
      */
     for(int k = 0; k < nopt; k++){
       c_dist2 = std::pow(des_x - x_Options[k], 2) +
         std::pow(des_y - y_Options[k], 2);
       c_dist = std::sqrt(c_dist2);
       distance_toDes[k] = c_dist;
     }
     // now we need to normalise the distances to something compatible with
     // the weighting for sample choice (xi – min(x)) / (max(x) – min(x))
     // find MIN
     double dist_min = distance_toDes[0];
     dist_min = cpp_min(distance_toDes);
     // find MAX
     double dist_max = distance_toDes[0];
     dist_max = cpp_max(distance_toDes);

     // using the find max and min found above we normalise all the distances
     // from possible choices in relation to the final destination
     for(int m = 0; m < nopt; m++){
       // first we have to invert it so ones closer to the destination are
       // preferred. ie larger numbers and greater weighting in the sample
       // function
       weights_toDes[m] = 1 - ((distance_toDes[m] - dist_min) /
         (dist_max - dist_min));
       // then combine them with the movement Matrix values
       // we can deal with some balancing issues here

       // only apply distance/point attraction to non-exploratory behaviours
       if(behave_Locations[i] == 1){
         move_Options[m] = move_Options[m];

       } else {
         // and use the different ways of balancing the influence of the destination
         switch(destinationTrans){
         case 0:
           move_Options[m] = move_Options[m] + destinationMod * weights_toDes[m];
           break;
         case 1:
           move_Options[m] = move_Options[m] + destinationMod * std::sqrt(weights_toDes[m]);
           break;
         case 2:
           move_Options[m] = move_Options[m] + destinationMod * std::pow(weights_toDes[m], 2);
           break;
         } // switch end
       } // if end
     }// for m end



     // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     // DEALING WITH AVOIDANCE POINTS
     // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     // current distances from all avoidance points
     cumulative_dist = 0;
     for(int mO = 0; mO < nopt; mO++){
       for(int avp = 0; avp < navp; avp++){

         c_dist2 = std::pow(avoidPoints_x[avp] - x_Options[mO], 2) +
           std::pow(avoidPoints_y[avp] - y_Options[mO], 2);
         c_dist = std::sqrt(c_dist2);
         cumulative_dist = cumulative_dist + c_dist;

       }
       distance_toAvoid[mO] = cumulative_dist;
     }

     // find MIN
     double distAvoid_min = distance_toAvoid[0];
     distAvoid_min = cpp_min(distance_toAvoid);
     // find MAX
     double distAvoid_max = distance_toAvoid[0];
     distAvoid_max = cpp_max(distance_toAvoid);

     // use the min and max cumdists to the avoidance points to modify the
     // move_Options draw
     for(int m = 0; m < nopt; m++){
       // This time we are not inverting it so ones farther from the avoidance
       // points are preferred.
       weights_toDes[m] = (distance_toAvoid[m] - dist_min) /
         (distAvoid_max - distAvoid_min);
       // then combine them with the movement Matrix values

       // and use the different ways of balancing the influence of avoidance
       switch(avoidTrans){
       case 0:
         move_Options[m] = move_Options[m] + avoidMod * weights_toDes[m];
         break;
       case 1:
         move_Options[m] = move_Options[m] + avoidMod * std::sqrt(weights_toDes[m]);
         break;
       case 2:
         move_Options[m] = move_Options[m] + avoidMod * std::pow(weights_toDes[m], 2);
         break;
       } // switch end
     }



     // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     // NOW THE AGENT MAKES ITS CHOICE
     // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     chosen = cpp_sample_options(move_Options);
     chosen_Options[i] = chosen;

     // make sure to update the direction of travel each move
     last_angle = taOptions[chosen];
     if(last_angle > 180){
       last_angle = last_angle - 180;
     } else if(last_angle < - 180){
       last_angle = last_angle + 180;
     }

     x_Locations[i]  = x_Options[chosen];
     y_Locations[i]  = y_Options[chosen];
     sl_Locations[i] = slOptions[chosen];
     ta_Locations[i] = last_angle;
     step_Locations[i] = i;

   }
   if (dolog) {std::cout << "| 100% HOORAY! \\o/" << std::endl;}

   Rcpp::List INPUT_basic = Rcpp::List::create(
     Rcpp::Named("in_startx") = startx,
     Rcpp::Named("in_starty") = starty,
     Rcpp::Named("in_timesteps") = timesteps,
     Rcpp::Named("in_ndes") = ndes,
     Rcpp::Named("in_nopt") = nopt
   );

   Rcpp::List INPUT_destination = Rcpp::List::create(
     Rcpp::Named("in_shelter_locs_x") = shelter_locs_x,
     Rcpp::Named("in_shelter_locs_y") = shelter_locs_y,
     Rcpp::Named("in_sSiteSize") = sSiteSize,
     Rcpp::Named("in_avoidPoints_x") = avoidPoints_x,
     Rcpp::Named("in_avoidPoints_y") = avoidPoints_y,
     Rcpp::Named("in_k_desRange") = k_desRange,
     Rcpp::Named("in_s_desRange") = s_desRange,
     Rcpp::Named("in_mu_desDir") = mu_desDir,
     Rcpp::Named("in_k_desDir") = k_desDir,
     Rcpp::Named("in_destinationTrans") = destinationTrans,
     Rcpp::Named("in_destinationMod") = destinationMod,
     Rcpp::Named("in_avoidTrans") = avoidTrans,
     Rcpp::Named("in_avoidMod") = avoidMod
   );

   Rcpp::List INPUT_movement = Rcpp::List::create(
     Rcpp::Named("in_k_step") = k_step,
     Rcpp::Named("in_s_step") = s_step,
     Rcpp::Named("in_mu_angle") = mu_angle,
     Rcpp::Named("in_k_angle") = k_angle,
     Rcpp::Named("in_rescale") = rescale,
     Rcpp::Named("in_b0_Options") = b0_Options,
     Rcpp::Named("in_b1_Options") = b1_Options,
     Rcpp::Named("in_b2_Options") = b2_Options
   );

   Rcpp::List INPUT_cycle = Rcpp::List::create(
     Rcpp::Named("in_rest_Cycle_A") = rest_Cycle_A,
     Rcpp::Named("in_rest_Cycle_M") = rest_Cycle_M,
     Rcpp::Named("in_rest_Cycle_PHI") = rest_Cycle_PHI,
     Rcpp::Named("in_rest_Cycle_TAU") = rest_Cycle_TAU,
     Rcpp::Named("in_addCycles") = addCycles,
     Rcpp::Named("in_add_Cycle_A") = add_Cycle_A,
     Rcpp::Named("in_add_Cycle_M") = add_Cycle_M,
     Rcpp::Named("in_add_Cycle_PHI") = add_Cycle_PHI,
     Rcpp::Named("in_add_Cycle_TAU") = add_Cycle_TAU
   );


   Rcpp::List INPUT_layers = Rcpp::List::create(
     Rcpp::Named("in_shelterMatrix") = shelterMatrix,
     Rcpp::Named("in_forageMatrix") = forageMatrix,
     Rcpp::Named("in_moveMatrix") = moveMatrix
   );

   Rcpp::List OUTPUT = Rcpp::List::create(
     // output the location data
     Rcpp::Named("loc_x") = x_Locations,
     Rcpp::Named("loc_y") = y_Locations,
     Rcpp::Named("loc_sl") = sl_Locations,
     Rcpp::Named("loc_ta") = ta_Locations,
     Rcpp::Named("loc_step") = step_Locations,
     Rcpp::Named("loc_step_rescale") = rescale,
     Rcpp::Named("loc_behave") = behave_Locations,
     // output for the chosen options at each step
     Rcpp::Named("loc_chosen") = chosen_Options,
     Rcpp::Named("loc_x_destinations") = des_x_Locations,
     Rcpp::Named("loc_y_destinations") = des_y_Locations,
     Rcpp::Named("loc_chosen_destinations") = des_chosen,
     // output for all the optionsALL
     Rcpp::Named("oall_x") = x_OptionsAll,
     Rcpp::Named("oall_y") = y_OptionsAll,
     Rcpp::Named("oall_step") = step_OptionsAll,
     Rcpp::Named("oall_stepLengths") = stepAll,
     // output for all dynamic desintations and the chosen shelters
     Rcpp::Named("destinations_list") = Rcpp::List::create(
       Rcpp::Named("des_desOptsx") = x_DesOptionsAll,
       Rcpp::Named("des_desOptsy") = y_DesOptionsAll,
       Rcpp::Named("des_desBehave") = behave_DesOptionsAll,
       Rcpp::Named("des_desOptsstep") = step_DesOptionsAll
     ),
     // output for the last options just to check
     // Rcpp::Named("opt_x_forOpts") = x_forageOptions,
     // Rcpp::Named("opt_y_forOpts") = y_forageOptions,
     // Rcpp::Named("opt_des_forOpts") = des_forageOptions,
     // Rcpp::Named("opt_chosen_forOpts") = chosenDes,
     // bring in all the input lists together
     Rcpp::Named("inputs_list") = Rcpp::List::create(
       Rcpp::Named("inputs_basic") = INPUT_basic,
       Rcpp::Named("inputs_destination") = INPUT_destination,
       Rcpp::Named("inputs_movement") = INPUT_movement,
       Rcpp::Named("inputs_cycle") = INPUT_cycle,
       Rcpp::Named("inputs_layers") = INPUT_layers)
     // Rcpp::Named("ol_x") = x_Options,
     // Rcpp::Named("ol_y") = y_Options,
     // Rcpp::Named("ol_moveVal1") = move_Options, // included to check probs used
     // Rcpp::Named("ol_step") = step_Options, // included to check is choice vector is the source of issues
     // Rcpp::Named("ol_c_dist2") = c_dist2,
     // Rcpp::Named("ol_c_dist") = std::sqrt(c_dist2),
     // Rcpp::Named("ol_dist2Des") = distance_toDes,
     // Rcpp::Named("ol_dist2DesInvert") = distInvert,
     // Rcpp::Named("ol_distWeights") = weights_toDes
   );
   return OUTPUT;

 }
