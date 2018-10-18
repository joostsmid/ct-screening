# writes the .cpp code utilized to do the calculations for the differentials
generate.cpp.code <- function(parms, filename = "STIreduced.cpp"){
  
  write(c("#include <stdio.h>","#include <stdlib.h>", "#include <math.h>", 
          "#define MATHLIB_STANDALONE 1", "#include <Rmath.h>", 
          "#include <vector>", " ", "using namespace std;", 
          "typedef vector<int> ivector;", "typedef vector<double> dvector;", " "), filename,  sep = "\n")
  
  #### beginning of the function
  write("extern \"C\" {", filename, append = TRUE, sep = "\n")
  
  ### declaration of the function
  write(c(paste("   void STImod(double *args, double *res)", sep = ""), "   {"), filename, append = TRUE, sep = "\n")
  write("  ", filename, append = TRUE, sep = "\n")
  
  
  #########################################################################################################################
  
  ### define auxiliary variables
  write("   int arg_counter = 0;", filename, append = TRUE, sep = "\n")
  write("   bool verbose = false;", filename, append = TRUE, sep = "\n") 
  write("  ", filename, append = TRUE, sep = "\n")
  
  
  #########################################################################################################################
  ### start to fill up parameters 
  
  ### time #JS150316
  write("   /// read in time", filename, append = TRUE, sep = "\n") #JS150316
  write("   double t = args[arg_counter];", filename, append = TRUE, sep = "\n") #JS150316
  write("   arg_counter++;", filename, append = TRUE, sep = "\n") #JS150316
  
  ### sex
  write(paste("   int sex_length = ", length(parms$sex),";", sep = ""), filename, append = TRUE, sep = "\n")
  write("  ", filename, append = TRUE, sep = "\n")
  
  
  ### nmb of activity classes
  write(paste("   int nmb_act_classes = ", parms$nJ ,";", sep = ""), filename, append = TRUE, sep = "\n")
  write("  ", filename, append = TRUE, sep = "\n")
  
  ### nmb of screening years #JS170316
  n_years_screening <- unname(dim(parms$chi_all)[3])
  write(paste("   int nmb_screen_years = ", n_years_screening ,";", sep = ""), filename, append = TRUE, sep = "\n") #JS170316
  write("  ", filename, append = TRUE, sep = "\n") #JS170316
  
  ### nmb of age classes
  write(paste("   int nmb_age_classes = ", parms$n.age.classes ,";", sep = ""), filename, append = TRUE, sep = "\n")
  write("  ", filename, append = TRUE, sep = "\n")
  
  ### aging rate
  write("   double aging_rate[nmb_age_classes];", filename, append = TRUE, sep = "\n")
  write("  ", filename, append = TRUE, sep = "\n")
  ### fill in 
  for(a in 1:parms$n.age.classes){
    write(paste("   aging_rate[", a -1  ,"] = ", parms$ar[a],";", sep = ""), filename, append = TRUE, sep = "\n")
  }
  write("  ", filename, append = TRUE, sep = "\n")
  
  ### virginity loss rate
  write("   double virg_loss[sex_length][nmb_age_classes];", filename, append = TRUE, sep = "\n")
  write("  ", filename, append = TRUE, sep = "\n")
  ### fill in 
  for(g in 1:length(parms$sex)){
    for(a in 1:parms$n.age.classes){
      write(paste("   virg_loss[", g -1  ,"][",a-1  ,"] = ", parms$vl[g,a],";", sep = ""), filename, append = TRUE, sep = "\n")
    }
  }
  write("  ", filename, append = TRUE, sep = "\n")
  
  ### rate at which virginity is lost to different sexual activity classes
  write("   double fold[sex_length][nmb_act_classes];", filename, append = TRUE, sep = "\n")
  write("  ", filename, append = TRUE, sep = "\n")
  
  ### fill this in 
  for(g in 1:length(parms$sex)){
    for(j in 1:parms$nJ){
      write(paste("   fold[", g -1  ,"][", j - 1  ,"] = ", parms$fold[g,j],";", sep = ""), filename, append = TRUE, sep = "\n")
    }
  }
  write("  ", filename, append = TRUE, sep = "\n")
  
  
  ### transmission rate (NOT PRE DEFINED, BY PASSED BY R)
  write("   double beta[sex_length];", filename, append = TRUE, sep = "\n") #JS010516
  write("   beta[0] = args[arg_counter];", filename, append = TRUE, sep = "\n") #JS010516
  ### add verbose code part
#   write("   if (verbose) {", filename, append = TRUE, sep = "\n")
#   write("      printf(\"transmission rate M per susceptible partnership is %f; \", beta[1]);", filename, append = TRUE, sep = "\n")
#   write("      printf(\"Arg count is at %d; \", arg_counter);", filename, append = TRUE, sep = "\n")
#   write("   }", filename, append = TRUE, sep = "\n")
  write("  ", filename, append = TRUE, sep = "\n")
  write("   arg_counter++;", filename, append = TRUE, sep = "\n") #JS010516
  
  write("   beta[1] = args[arg_counter];", filename, append = TRUE, sep = "\n") #JS010516
#   write("   if (verbose) {", filename, append = TRUE, sep = "\n")
#   write("      printf(\"transmission rate F per susceptible partnership is %f; \", beta[2]);", filename, append = TRUE, sep = "\n")
#   write("      printf(\"Arg count is at %d; \", arg_counter);", filename, append = TRUE, sep = "\n")
#   write("   }", filename, append = TRUE, sep = "\n")
  write("  ", filename, append = TRUE, sep = "\n")
  write("   arg_counter++;", filename, append = TRUE, sep = "\n") #JS010516
  write("  ", filename, append = TRUE, sep = "\n") #JS010516
  
  ### contact rate 
  write("   double ct[sex_length][nmb_act_classes][nmb_age_classes];", filename, append = TRUE, sep = "\n")
  write("  ", filename, append = TRUE, sep = "\n")
  
  ### fill in contact rates
  for(g in 1:length(parms$sex)){
    for(j in 1:parms$nJ){
      for(a in 1:parms$n.age.classes){
        write(paste("   ct[", g -1  ,"][",j -1 ,"][",a-1,"] = ", parms$ct[g,j,a],";", sep = ""), filename, append = TRUE, sep = "\n") #JS080416
      }
    }
  }
  write("  ", filename, append = TRUE, sep = "\n")
  
  ### age mixing matrix
  write("   double rho_age[sex_length][nmb_age_classes][nmb_age_classes];", filename, append = TRUE, sep = "\n")
  write("  ", filename, append = TRUE, sep = "\n")
  
  ### fill in the matrix
  for(g in 1:length(parms$sex)){
    for(a in 1:parms$n.age.classes){
      for(apr in 1:parms$n.age.classes){
        write(paste("   rho_age[", g -1  ,"][", a -1 ,"][",apr-1,"] = ", parms$rho_age[g,a,apr],";", sep = ""), filename, append = TRUE, sep = "\n")
      }
    }
  }
  write("  ", filename, append = TRUE, sep = "\n")
  
  ### balancing exponent theta
  write("   double theta[sex_length];", filename, append = TRUE, sep = "\n")
  write("  ", filename, append = TRUE, sep = "\n")
  
  ### fill in 
  for(g in 1:length(parms$sex)){
    write(paste("   theta[", g -1  ,"] = ", parms$theta[g],";", sep = ""), filename, append = TRUE, sep = "\n")
  }
  write("  ", filename, append = TRUE, sep = "\n")
  
  ### transmission rate (NOT PRE DEFINED, BY PASSED BY R)
  write("   double epsilon = args[arg_counter];", filename, append = TRUE, sep = "\n")
  write("   arg_counter++;", filename, append = TRUE, sep = "\n")
  write("  ", filename, append = TRUE, sep = "\n")
  
#   ### add verbose code part
#   write("   if (verbose) {", filename, append = TRUE, sep = "\n")
#   write("      printf(\"have implemented up to epsilon rate\");", filename, append = TRUE, sep = "\n")
#   write("      printf(\"value of epsilon is %f; \", epsilon);", filename, append = TRUE, sep = "\n")
#   write("      printf(\"Arg count is at %d; \", arg_counter);", filename, append = TRUE, sep = "\n")
#   write("   }", filename, append = TRUE, sep = "\n")
#   write("  ", filename, append = TRUE, sep = "\n")
  
  ### switchting rates between activity classes, depending on gender and age
  write("   double sw[sex_length][nmb_act_classes][nmb_act_classes][nmb_age_classes];", filename, append = TRUE, sep = "\n")
  write("     ", filename, append = TRUE, sep = "\n")
  
  ### fill in switching rates
  for (g in 1:length(parms$sex)) {
    for (j in 1:parms$nJ) {
      for (jpr in 1:parms$nJ) {
        for (a in 1:parms$n.age.classes) {
          write(paste("   sw[", g -1  ,"][", j -1 ,"][", jpr -1 ,"][",a-1,"] = ", parms$sw[g,j,jpr,a],";", sep = ""), filename, append = TRUE, sep = "\n")
        }
      }
    }
  }
  
  ### clearance rate of infection
  write("   double gamma[sex_length];", filename, append = TRUE, sep = "\n") #JS010516
  write("   gamma[0] = args[arg_counter];", filename, append = TRUE, sep = "\n") #JS010516
#   ### add verbose code part
#   write("   if (verbose) {", filename, append = TRUE, sep = "\n")
#   write("      printf(\"clearance rate M per susceptible partnership is %f; \", gamma[0]);", filename, append = TRUE, sep = "\n")
#   write("      printf(\"Arg count is at %d; \", arg_counter);", filename, append = TRUE, sep = "\n")
#   write("   }", filename, append = TRUE, sep = "\n")
  write("  ", filename, append = TRUE, sep = "\n")
  write("   arg_counter++;", filename, append = TRUE, sep = "\n") #JS010516
  
  write("   gamma[1] = args[arg_counter];", filename, append = TRUE, sep = "\n") #JS010516
#   ### add verbose code part
#   write("   if (verbose) {", filename, append = TRUE, sep = "\n")
#   write("      printf(\"clearance rate F per susceptible partnership is %f; \", gamma[1]);", filename, append = TRUE, sep = "\n")
#   write("      printf(\"Arg count is at %d; \", arg_counter);", filename, append = TRUE, sep = "\n")
#   write("   }", filename, append = TRUE, sep = "\n")
  write("  ", filename, append = TRUE, sep = "\n")
  write("   arg_counter++;", filename, append = TRUE, sep = "\n") #JS010516
  write("  ", filename, append = TRUE, sep = "\n") #JS010516
  
  ### factor by which susceptibility to reinfection is reduced
  write("   double kappa[sex_length];", filename, append = TRUE, sep = "\n") #JS010516
  write("   kappa[0] = args[arg_counter];", filename, append = TRUE, sep = "\n") #JS010516
  write("   arg_counter++;", filename, append = TRUE, sep = "\n") #JS010516
  write("   kappa[1] = args[arg_counter];", filename, append = TRUE, sep = "\n") #JS010516
  write("   arg_counter++;", filename, append = TRUE, sep = "\n")
  write("  ", filename, append = TRUE, sep = "\n")
  
  ### efficacy of treatment
  write("   double omega[2];", filename, append = TRUE, sep = "\n") 
  write("   omega[0] = args[arg_counter];", filename, append = TRUE, sep = "\n")
  write("   arg_counter++;", filename, append = TRUE, sep = "\n") 
  write("   omega[1] = args[arg_counter];", filename, append = TRUE, sep = "\n") 
  write("   arg_counter++;", filename, append = TRUE, sep = "\n") 
  write("  ", filename, append = TRUE, sep = "\n")

  
  ### fraction symptomatic
  write("   double fsymp[sex_length];", filename, append = TRUE, sep = "\n") #JS010516
  write("   fsymp[0] = args[arg_counter];", filename, append = TRUE, sep = "\n") #JS010516
  write("   arg_counter++;", filename, append = TRUE, sep = "\n") #JS010516
  write("   fsymp[1] = args[arg_counter];", filename, append = TRUE, sep = "\n") #JS010516
  write("   arg_counter++;", filename, append = TRUE, sep = "\n")
  write("  ", filename, append = TRUE, sep = "\n")
  
  ### treatment rate
  write("   double treat;", filename, append = TRUE, sep = "\n") #JS250716
  write("   treat = args[arg_counter];", filename, append = TRUE, sep = "\n") #JS250716
  write("   arg_counter++;", filename, append = TRUE, sep = "\n") #JS250716
  write("  ", filename, append = TRUE, sep = "\n")
  
  ### screening uptake in susceptibles/recovereds wrt screening uptake infecteds/exposed, for different years
  write("   double eta[2];", filename, append = TRUE, sep = "\n")
  for(y in 1:2){
    write(paste("   eta[", y - 1  ,"] = args[arg_counter];", sep = ""), filename, append = TRUE, sep = "\n")
    write("   arg_counter++;", filename, append = TRUE, sep = "\n") 
  }
  write("  ", filename, append = TRUE, sep = "\n")
  
  ### screening uptake, NCSP data (total number of tests/year in total population)
  write("   double chi_all[sex_length][nmb_age_classes][nmb_screen_years];", filename, append = TRUE, sep = "\n") #JS170117
  write("  ", filename, append = TRUE, sep = "\n") #JS170117
  for(g in 1:length(parms$sex)){ #JS170117
    for(a in 1:parms$n.age.classes){ #JS170117
      for(y in 1:n_years_screening){ #JS170117
        write(paste("   chi_all[", g -1  ,"][",a-1,"][", y -1  ,"] = ", parms$chi_all[g,a,y],";", sep = ""), filename, append = TRUE, sep = "\n")
      } #JS170117
    } #JS170117
  } #JS170117
  write("  ", filename, append = TRUE, sep = "\n")
  
  ### assumed period of testing before intervention
  write(paste("   double testperiod_in_burnin = ", parms$testperiod_in_burnin,";", sep = ""), filename, append = TRUE, sep = "\n")
  write("  ", filename, append = TRUE, sep = "\n")
  
  ### burnintime
  write("   double burnintime;", filename, append = TRUE, sep = "\n")
  write(paste("   burnintime = ", parms$burnintime, ";", sep = ""), filename, append = TRUE, sep = "\n")
  write("  ", filename, append = TRUE, sep = "\n")
  
  ###########################################  DEFINITION OF VARIABLES ##############################################################################
  
  write("   ///////////////////////////////////////////////////////////////////////////////", filename, append = TRUE, sep = "\n")
  write("   /// define the essential variables ///////", filename, append = TRUE, sep = "\n")
  write("   ///////////////////////////////////////////////////////////////////////////////", filename, append = TRUE, sep = "\n")
  write("  ", filename, append = TRUE, sep = "\n")
  
  write("   double U[sex_length][nmb_age_classes]; // virgins", filename, append = TRUE, sep = "\n")
  write("   double S[sex_length][nmb_act_classes][nmb_age_classes]; // susceptibles", filename, append = TRUE, sep = "\n")
  write("   double I_A[sex_length][nmb_act_classes][nmb_age_classes]; // asymptomatic infecteds", filename, append = TRUE, sep = "\n")
  write("   double I_S[sex_length][nmb_act_classes][nmb_age_classes]; // symptomatic infecteds", filename, append = TRUE, sep = "\n")
  write("   double R[sex_length][nmb_act_classes][nmb_age_classes]; // recovereds", filename, append = TRUE, sep = "\n")
  write("   double I_A2[sex_length][nmb_act_classes][nmb_age_classes]; // asymptomatic infecteds", filename, append = TRUE, sep = "\n")
  write("   double I_S2[sex_length][nmb_act_classes][nmb_age_classes]; // symptomatic infecteds", filename, append = TRUE, sep = "\n")
  write("   double N[sex_length][nmb_act_classes][nmb_age_classes]; // sum of all", filename, append = TRUE, sep = "\n")
  write("  ", filename, append = TRUE, sep = "\n")
  
  write("   int tensor_entry_numbers = sex_length * nmb_act_classes * nmb_age_classes;", filename, append = TRUE, sep = "\n")
  write("  ", filename, append = TRUE, sep = "\n")
  
  write("   /// define new counter ", filename, append = TRUE, sep = "\n")
  write("   int counter = arg_counter;", filename, append = TRUE, sep = "\n")
  write("  ", filename, append = TRUE, sep = "\n")
  
  ### read in virgins
  write("   /// read in U: virgins", filename, append = TRUE, sep = "\n")
  write("   for(int g = 0; g < sex_length; g++){", filename, append = TRUE, sep = "\n")
  write("      for (int a = 0; a < nmb_age_classes; a++) {", filename, append = TRUE, sep = "\n")
  write("         U [g][a] = args[counter];", filename, append = TRUE, sep = "\n")
  write("      ", filename, append = TRUE, sep = "\n")
  write("          if(U [g][a] < 0){", filename, append = TRUE, sep = "\n")
  write("            if(verbose){", filename, append = TRUE, sep = "\n")
  write("               printf(\"NEGATIVE U entry! U[%d][%d]: %f; \", g, a, U[g][a]);", filename, append = TRUE, sep = "\n")
  write("            }", filename, append = TRUE, sep = "\n")
  write("         }", filename, append = TRUE, sep = "\n")
  write("      ", filename, append = TRUE, sep = "\n")
  write("         counter++;", filename, append = TRUE, sep = "\n")
  write("      }", filename, append = TRUE, sep = "\n")
  write("   } ", filename, append = TRUE, sep = "\n")
  write("  ", filename, append = TRUE, sep = "\n")
  
  ### read in initial conditions
  write("   /// read in the initial conditions passed on by the user call", filename, append = TRUE, sep = "\n")
  write("   for(int g = 0; g < sex_length; g++){", filename, append = TRUE, sep = "\n")
  write("      for (int j = 0; j < nmb_act_classes; j++) {", filename, append = TRUE, sep = "\n")
  write("         for (int a = 0; a < nmb_age_classes; a++) {", filename, append = TRUE, sep = "\n")
  write("   ", filename, append = TRUE, sep = "\n")
  write("            /// proceed to fill the tensors", filename, append = TRUE, sep = "\n")
  
  write("            S [g][j][a] = args[counter];", filename, append = TRUE, sep = "\n")
  write("   ", filename, append = TRUE, sep = "\n")
  write("            /// check reasonableness for each tensor", filename, append = TRUE, sep = "\n")
  write("            if(S [g][j][a] < 0){", filename, append = TRUE, sep = "\n")
  write("               if(verbose){", filename, append = TRUE, sep = "\n")
  write("                  printf(\"NEGATIVE S entry! S[%d][%d][%d]: %f; \", g, j, a, S[g][j][a]);", filename, append = TRUE, sep = "\n")
  write("               }", filename, append = TRUE, sep = "\n")
  write("               S [g][j][a] = 0;", filename, append = TRUE, sep = "\n")
  write("            }", filename, append = TRUE, sep = "\n")
  write("            ", filename, append = TRUE, sep = "\n")
  
  write("            I_A [g][j][a] = args[counter + 1*tensor_entry_numbers];", filename, append = TRUE, sep = "\n")
  write("            ", filename, append = TRUE, sep = "\n")
  write("            /// check reasonableness for each tensor", filename, append = TRUE, sep = "\n")
  write("            if(I_A [g][j][a] < 0){", filename, append = TRUE, sep = "\n")
  write("               if(verbose){", filename, append = TRUE, sep = "\n")
  write("                  printf(\"NEGATIVE I_A entry! I_A[%d][%d][%d]: %f; \", g, j, a, I_A[g][j][a]);", filename, append = TRUE, sep = "\n")
  write("               }", filename, append = TRUE, sep = "\n")
  write("               ", filename, append = TRUE, sep = "\n")
  write("               I_A [g][j][a] = 0;", filename, append = TRUE, sep = "\n")
  write("            }", filename, append = TRUE, sep = "\n")
  write("   ", filename, append = TRUE, sep = "\n")
  
  write("            I_S [g][j][a] = args[counter + 2*tensor_entry_numbers];", filename, append = TRUE, sep = "\n")
  write("            ", filename, append = TRUE, sep = "\n")
  write("            /// check reasonableness for each tensor", filename, append = TRUE, sep = "\n")
  write("            if(I_S [g][j][a] < 0){", filename, append = TRUE, sep = "\n")
  write("               if(verbose){", filename, append = TRUE, sep = "\n")
  write("                  printf(\"NEGATIVE I_S entry! I_S[%d][%d][%d]: %f; \", g, j, a, I_S[g][j][a]);", filename, append = TRUE, sep = "\n")
  write("               }", filename, append = TRUE, sep = "\n")
  write("               ", filename, append = TRUE, sep = "\n")
  write("               I_S [g][j][a] = 0;", filename, append = TRUE, sep = "\n")
  write("            }", filename, append = TRUE, sep = "\n")
  write("   ", filename, append = TRUE, sep = "\n")
  
  write("            R [g][j][a] = args[counter + 3*tensor_entry_numbers];", filename, append = TRUE, sep = "\n")
  write("    ", filename, append = TRUE, sep = "\n")
  write("            /// check reasonableness for each tensor", filename, append = TRUE, sep = "\n")
  write("            if(R [g][j][a] < 0){", filename, append = TRUE, sep = "\n")
  write("               if(verbose){", filename, append = TRUE, sep = "\n")
  write("                  printf(\"NEGATIVE R entry! S[%d][%d][%d]: %f; \", g, j, a, R[g][j][a]);", filename, append = TRUE, sep = "\n")
  write("               }", filename, append = TRUE, sep = "\n")
  write("               ", filename, append = TRUE, sep = "\n")
  write("               R [g][j][a] = 0;", filename, append = TRUE, sep = "\n")
  write("            }", filename, append = TRUE, sep = "\n")
  write("   ", filename, append = TRUE, sep = "\n")
  
  write("            I_A2 [g][j][a] = args[counter + 4*tensor_entry_numbers];", filename, append = TRUE, sep = "\n")
  write("            ", filename, append = TRUE, sep = "\n")
  write("            /// check reasonableness for each tensor", filename, append = TRUE, sep = "\n")
  write("            if(I_A2 [g][j][a] < 0){", filename, append = TRUE, sep = "\n")
  write("               if(verbose){", filename, append = TRUE, sep = "\n")
  write("                  printf(\"NEGATIVE I_A2 entry! I_A2[%d][%d][%d]: %f; \", g, j, a, I_A2[g][j][a]);", filename, append = TRUE, sep = "\n")
  write("               }", filename, append = TRUE, sep = "\n")
  write("               ", filename, append = TRUE, sep = "\n")
  write("               I_A2 [g][j][a] = 0;", filename, append = TRUE, sep = "\n")
  write("            }", filename, append = TRUE, sep = "\n")
  write("   ", filename, append = TRUE, sep = "\n")
  
  write("            I_S2 [g][j][a] = args[counter + 5*tensor_entry_numbers];", filename, append = TRUE, sep = "\n")
  write("            ", filename, append = TRUE, sep = "\n")
  write("            /// check reasonableness for each tensor", filename, append = TRUE, sep = "\n")
  write("            if(I_S2 [g][j][a] < 0){", filename, append = TRUE, sep = "\n")
  write("               if(verbose){", filename, append = TRUE, sep = "\n")
  write("                  printf(\"NEGATIVE I_S2 entry! I_S2[%d][%d][%d]: %f; \", g, j, a, I_S2[g][j][a]);", filename, append = TRUE, sep = "\n")
  write("               }", filename, append = TRUE, sep = "\n")
  write("               ", filename, append = TRUE, sep = "\n")
  write("               I_S2 [g][j][a] = 0;", filename, append = TRUE, sep = "\n")
  write("            }", filename, append = TRUE, sep = "\n")
  write("   ", filename, append = TRUE, sep = "\n")
  
  write("            N [g][j][a] = S [g][j][a] + I_A [g][j][a] + I_S [g][j][a] + R [g][j][a] + I_A2 [g][j][a] + I_S2 [g][j][a];", filename, append = TRUE, sep = "\n")
  write("   ", filename, append = TRUE, sep = "\n")
  write("            if(N [g][j][a] < 0){", filename, append = TRUE, sep = "\n")
  write("    ", filename, append = TRUE, sep = "\n")
  write("               if(verbose){", filename, append = TRUE, sep = "\n")
  write("                  printf(\"NEGATIVE N entry! N[%d][%d][%d]: %f; \", g, j, a, N[g][j][a]);", filename, append = TRUE, sep = "\n")
  write("               }", filename, append = TRUE, sep = "\n")
  write("               ", filename, append = TRUE, sep = "\n")
  write("               N [g][j][a] = 0;", filename, append = TRUE, sep = "\n")
  write("            }", filename, append = TRUE, sep = "\n")
  write("    ", filename, append = TRUE, sep = "\n")
#   write("            if(verbose){", filename, append = TRUE, sep = "\n")
#   write("               printf(\"S[%d][%d][%d]: %f; \", g, j, a, S[g][j][a]);", filename, append = TRUE, sep = "\n")
#   write("               printf(\"I_A[%d][%d][%d]: %f; \", g, j, a, I_A[g][j][a]);", filename, append = TRUE, sep = "\n")
#   write("               printf(\"I_S[%d][%d][%d]: %f; \", g, j, a, I_S[g][j][a]);", filename, append = TRUE, sep = "\n")
#   write("               printf(\"R[%d][%d][%d]: %f; \", g, j, a, R[g][j][a]);", filename, append = TRUE, sep = "\n")
#   write("               printf(\"I_A2[%d][%d][%d]: %f; \", g, j, a, I_A2[g][j][a]);", filename, append = TRUE, sep = "\n")
#   write("               printf(\"I_S2[%d][%d][%d]: %f; \", g, j, a, I_S2[g][j][a]);", filename, append = TRUE, sep = "\n")
#   write("            }", filename, append = TRUE, sep = "\n")
  write("               ", filename, append = TRUE, sep = "\n")
  write("            counter++;", filename, append = TRUE, sep = "\n")
  write("         }", filename, append = TRUE, sep = "\n")
  write("      }", filename, append = TRUE, sep = "\n")
  write("   } ", filename, append = TRUE, sep = "\n")
  
  ###########################################  PERFORM CALCULATIONS ##############################################################################
  
  write("  ", filename, append = TRUE, sep = "\n")
  write("   ///////////////////////////////////////////////////////////////////////////////", filename, append = TRUE, sep = "\n")
  write("   /// test data", filename, append = TRUE, sep = "\n")
  write("   ///////////////////////////////////////////////////////////////////////////////", filename, append = TRUE, sep = "\n")
  write("    ", filename, append = TRUE, sep = "\n")

  write("   /// total number of screen tests in compartments the entire population (time dependent)", filename, append = TRUE, sep = "\n")
  write("   double screen_all[sex_length][nmb_age_classes];", filename, append = TRUE, sep = "\n")
  write("   ", filename, append = TRUE, sep = "\n")
  write("   for(int g = 0; g < sex_length; g++){  ", filename, append = TRUE, sep = "\n")
  write("     for (int a = 0; a < nmb_age_classes; a++) {", filename, append = TRUE, sep = "\n")
  write("       if(t < burnintime-testperiod_in_burnin){", filename, append = TRUE, sep = "\n")
  write("         screen_all[g][a]= 0;", filename, append = TRUE, sep = "\n")
  write("       } else if (t < burnintime) {", filename, append = TRUE, sep = "\n")
  write("         double I_S_ga = 0;", filename, append = TRUE, sep = "\n")
  write("         for (int j = 0; j < nmb_act_classes; j++) {", filename, append = TRUE, sep = "\n")
  write("           I_S_ga = I_S_ga + I_S[g][j][a]+ I_S2[g][j][a];", filename, append = TRUE, sep = "\n")
  write("         }", filename, append = TRUE, sep = "\n")
  write("           if (treat * I_S_ga < chi_all[g][a][0]){", filename, append = TRUE, sep = "\n")
  write("              screen_all[g][a]= (t-(burnintime - testperiod_in_burnin)) / testperiod_in_burnin * (chi_all[g][a][0] - treat * I_S_ga);", filename, append = TRUE, sep = "\n")
  write("           } else {", filename, append = TRUE, sep = "\n")
  write("               screen_all[g][a]= 0;", filename, append = TRUE, sep = "\n")
  write("           }", filename, append = TRUE, sep = "\n")
  write("      } else if (t < (double) nmb_screen_years + burnintime) {", filename, append = TRUE, sep = "\n")
  write("        double I_S_ga = 0;", filename, append = TRUE, sep = "\n")
  write("        for (int j = 0; j < nmb_act_classes; j++) {", filename, append = TRUE, sep = "\n")
  write("          I_S_ga = I_S_ga + I_S[g][j][a] + I_S2[g][j][a];", filename, append = TRUE, sep = "\n")
  write("        }", filename, append = TRUE, sep = "\n")
  write("           if (treat * I_S_ga < chi_all[g][a][(int) floor(t - burnintime)]){", filename, append = TRUE, sep = "\n")
  write("              screen_all[g][a]= chi_all[g][a][(int) floor(t - burnintime)] - treat * I_S_ga;", filename, append = TRUE, sep = "\n")
  write("           } else {", filename, append = TRUE, sep = "\n")
  write("               screen_all[g][a]= 0;", filename, append = TRUE, sep = "\n")
  write("           }", filename, append = TRUE, sep = "\n")
  write("       } else {", filename, append = TRUE, sep = "\n")
  write("         double I_S_ga = 0;", filename, append = TRUE, sep = "\n")
  write("         for (int j = 0; j < nmb_act_classes; j++) {", filename, append = TRUE, sep = "\n")
  write("           I_S_ga = I_S_ga + I_S[g][j][a] + I_S2[g][j][a];", filename, append = TRUE, sep = "\n")
  write("         }", filename, append = TRUE, sep = "\n")
  write("           if (treat * I_S_ga < chi_all[g][a][nmb_screen_years-1]){", filename, append = TRUE, sep = "\n")
  write("              screen_all[g][a]= chi_all[g][a][nmb_screen_years-1]  - treat * I_S_ga;", filename, append = TRUE, sep = "\n")
  write("           } else {", filename, append = TRUE, sep = "\n")
  write("               screen_all[g][a]= 0;", filename, append = TRUE, sep = "\n")
  write("           }", filename, append = TRUE, sep = "\n")
  write("       }", filename, append = TRUE, sep = "\n")
  write("     }", filename, append = TRUE, sep = "\n")
  write("   }", filename, append = TRUE, sep = "\n")
  write("   ", filename, append = TRUE, sep = "\n")
  
  write("     /// eta_exp as exponential relation of number of offered tests", filename, append = TRUE, sep = "\n")
  write("   double eta_exp[sex_length][nmb_age_classes];", filename, append = TRUE, sep = "\n")
  write("   for(int g = 0; g < sex_length; g++){  ", filename, append = TRUE, sep = "\n")
  write("     for (int a = 0; a < nmb_age_classes; a++) {", filename, append = TRUE, sep = "\n")
  write("         eta_exp[g][a] =1 + (eta[0]-1) * exp(-eta[1]*screen_all[g][a]);", filename, append = TRUE, sep = "\n")
  # write("         eta_exp[g][a] =1 + eta[0] - exp(eta[1]*screen_all[g][a]);", filename, append = TRUE, sep = "\n") # JS231017: eta convex assumed
  # write("         if (eta_exp[g][a] < 1){", filename, append = TRUE, sep = "\n") # JS231017: eta convex assumed
  # write("              eta_exp[g][a] = 1;", filename, append = TRUE, sep = "\n")
  # write("         }", filename, append = TRUE, sep = "\n")
  write("     }", filename, append = TRUE, sep = "\n")
  write("   }", filename, append = TRUE, sep = "\n")
  
  write("   /// total number of screen tests for infected persons per person eligible for screening (time dependent, heterogeneity between suscpetibles/infecteds accounted for)", filename, append = TRUE, sep = "\n")
  write("   double screen[sex_length][nmb_age_classes];", filename, append = TRUE, sep = "\n")
  write("   for(int g = 0; g < sex_length; g++){", filename, append = TRUE, sep = "\n")
  write("       for (int a = 0; a < nmb_age_classes; a++) {", filename, append = TRUE, sep = "\n")
  write("         double denom = 0;", filename, append = TRUE, sep = "\n")
  write("         for (int j = 0; j < nmb_act_classes; j++) {", filename, append = TRUE, sep = "\n")
  write("             denom = denom + S[g][j][a] + R[g][j][a] + eta_exp[g][a] * (I_A[g][j][a] + I_A2[g][j][a]);", filename, append = TRUE, sep = "\n")
  write("         }", filename, append = TRUE, sep = "\n")
  write("         screen[g][a] = eta_exp[g][a] * screen_all[g][a] / denom;", filename, append = TRUE, sep = "\n")
  write("       }", filename, append = TRUE, sep = "\n")
  write("   }", filename, append = TRUE, sep = "\n")
  
  write("     ///////////////////////////////////////////////////////////////////////////////", filename, append = TRUE, sep = "\n")
  write("     /// define the differentials (or rates) fo variables", filename, append = TRUE, sep = "\n")
  write("   ///////////////////////////////////////////////////////////////////////////////", filename, append = TRUE, sep = "\n")

  write("   double dU[sex_length][nmb_age_classes]; // virgins", filename, append = TRUE, sep = "\n")
  write("   double dS[sex_length][nmb_act_classes][nmb_age_classes]; // susceptibles", filename, append = TRUE, sep = "\n")
  write("   double dI_A[sex_length][nmb_act_classes][nmb_age_classes]; // asymp infecteds", filename, append = TRUE, sep = "\n")
  write("   double dI_S[sex_length][nmb_act_classes][nmb_age_classes]; // symp infecteds", filename, append = TRUE, sep = "\n")
  write("   double dR[sex_length][nmb_act_classes][nmb_age_classes]; // recovereds", filename, append = TRUE, sep = "\n")
  write("   double dI_A2[sex_length][nmb_act_classes][nmb_age_classes]; // asymp infecteds", filename, append = TRUE, sep = "\n")
  write("   double dI_S2[sex_length][nmb_act_classes][nmb_age_classes]; // symp infecteds", filename, append = TRUE, sep = "\n")
  write("   double dD[sex_length][nmb_act_classes][nmb_age_classes]; // diagnosed", filename, append = TRUE, sep = "\n")
  write("   double dInc[sex_length][nmb_act_classes][nmb_age_classes]; // incidence", filename, append = TRUE, sep = "\n")

  write("    ", filename, append = TRUE, sep = "\n")
  write("   /// calculate rho /////", filename, append = TRUE, sep = "\n")
  write("   double rho[sex_length][nmb_act_classes][nmb_act_classes][nmb_age_classes][nmb_age_classes];", filename, append = TRUE, sep = "\n")
  write("    ", filename, append = TRUE, sep = "\n")
  write("   /// auxiliary variables for rho", filename, append = TRUE, sep = "\n")
  write("   double kronecker = 0;", filename, append = TRUE, sep = "\n")
  write("   int gender_switch = 0;", filename, append = TRUE, sep = "\n")
  write("   double contacts_across_actclass = 0;", filename, append = TRUE, sep = "\n")
  write("    ", filename, append = TRUE, sep = "\n")
  write("   for(int g = 0; g < sex_length; g++){", filename, append = TRUE, sep = "\n")
  write("     for (int j = 0; j < nmb_act_classes; j++) {", filename, append = TRUE, sep = "\n")
  write("       for (int jpr = 0; jpr < nmb_act_classes; jpr++) {", filename, append = TRUE, sep = "\n")
  write("         for (int a = 0; a < nmb_age_classes; a++) {", filename, append = TRUE, sep = "\n")
  write("           for (int apr = 0; apr < nmb_age_classes; apr++) {", filename, append = TRUE, sep = "\n")
  write("             ", filename, append = TRUE, sep = "\n")
  write("             /// calc kronecker delta", filename, append = TRUE, sep = "\n")
  write("             j == jpr ? kronecker = 1 : kronecker = 0;", filename, append = TRUE, sep = "\n")
  write("             g == 0 ? gender_switch =1 : gender_switch = 0;", filename, append = TRUE, sep = "\n")
  write("             ", filename, append = TRUE, sep = "\n")
  write("             /// buld sum over all possible contacts", filename, append = TRUE, sep = "\n")
  write("             for (int l = 0; l < nmb_act_classes; l++) {", filename, append = TRUE, sep = "\n")
  write("               contacts_across_actclass += N[gender_switch][l][apr]*ct[gender_switch][l][apr];", filename, append = TRUE, sep = "\n")
  write("            }", filename, append = TRUE, sep = "\n")
  write("    ", filename, append = TRUE, sep = "\n")
  write("             /// set the rho value to zero if the contacts across are zero", filename, append = TRUE, sep = "\n")
  write("             if(contacts_across_actclass == 0){", filename, append = TRUE, sep = "\n")
  write("               ", filename, append = TRUE, sep = "\n")
  write("               rho [g][j][jpr][a][apr] = 0;", filename, append = TRUE, sep = "\n")
  write("               ", filename, append = TRUE, sep = "\n")
  write("             }", filename, append = TRUE, sep = "\n")
  write("             else{", filename, append = TRUE, sep = "\n")
  write("               /// fill in rho", filename, append = TRUE, sep = "\n")
  write("               rho [g][j][jpr][a][apr] = rho_age[g][a][apr]*(epsilon*kronecker + (1-epsilon)*(N[gender_switch][jpr][apr]*ct[gender_switch][jpr][apr])/contacts_across_actclass);", filename, append = TRUE, sep = "\n")
  write("               ", filename, append = TRUE, sep = "\n")
#   write("               if(verbose){", filename, append = TRUE, sep = "\n")
#   write("                 printf(\"rho_age[%d][%d][%d]: %f; \", g, a, apr, rho_age[g][a][apr]);", filename, append = TRUE, sep = "\n")
#   write("               }", filename, append = TRUE, sep = "\n")
  write("               ", filename, append = TRUE, sep = "\n")
  write("               if(rho [g][j][jpr][a][apr] < 0){", filename, append = TRUE, sep = "\n")
  write("                  if(verbose){", filename, append = TRUE, sep = "\n")
  write("                   printf(\"FOUND A NEGATIVE ENTRY IN RHO. rho [g][j][jpr][a][apr] is: %f \", rho [g][j][jpr][a][apr]);", filename, append = TRUE, sep = "\n")
  write("                   printf(\"N[%d][%d][%d]: %f; \", gender_switch,jpr, apr, N[gender_switch][jpr][apr]);", filename, append = TRUE, sep = "\n")
  write("                   printf(\"ct[%d][%d][%d]: %f; \", gender_switch,jpr, apr,ct[gender_switch][jpr][apr]);", filename, append = TRUE, sep = "\n")
  write("                   printf(\"rho_age[%d][%d][%d]: %f; \", g, j, a, rho_age[g][a][apr]);", filename, append = TRUE, sep = "\n")
  write("                   printf(\"contacts_across_actclass: %f \", contacts_across_actclass);", filename, append = TRUE, sep = "\n")
  write("                 }", filename, append = TRUE, sep = "\n")
  write("               }", filename, append = TRUE, sep = "\n")
  write("    ", filename, append = TRUE, sep = "\n")
  write("             }", filename, append = TRUE, sep = "\n")
  write("    ", filename, append = TRUE, sep = "\n")
  write("             contacts_across_actclass = 0.0;", filename, append = TRUE, sep = "\n")
  write("           }", filename, append = TRUE, sep = "\n")
  write("         }", filename, append = TRUE, sep = "\n")
  write("       }", filename, append = TRUE, sep = "\n")
  write("     }", filename, append = TRUE, sep = "\n")
  write("   } ", filename, append = TRUE, sep = "\n")
  
  write("    ", filename, append = TRUE, sep = "\n")
  write("   /// calculate contact rate ctc", filename, append = TRUE, sep = "\n")
  write("    ", filename, append = TRUE, sep = "\n")
  write("   double ctc[sex_length][nmb_act_classes][nmb_act_classes][nmb_age_classes][nmb_age_classes];", filename, append = TRUE, sep = "\n")
  write("   double cont_frac = 0;", filename, append = TRUE, sep = "\n")
  write("    ", filename, append = TRUE, sep = "\n")
  write("   for(int g = 0; g < sex_length; g++){", filename, append = TRUE, sep = "\n")
  write("     for (int j = 0; j < nmb_act_classes; j++) {", filename, append = TRUE, sep = "\n")
  write("       for (int jpr = 0; jpr < nmb_act_classes; jpr++) {", filename, append = TRUE, sep = "\n")
  write("         for (int a = 0; a < nmb_age_classes; a++) {", filename, append = TRUE, sep = "\n")
  write("           for (int apr = 0; apr < nmb_age_classes; apr++) {", filename, append = TRUE, sep = "\n")
  write("             ", filename, append = TRUE, sep = "\n")
  write("             /// gender switch", filename, append = TRUE, sep = "\n")
  write("             g == 0 ? gender_switch = 1 : gender_switch = 0;", filename, append = TRUE, sep = "\n")
  write("    ", filename, append = TRUE, sep = "\n")
  
  write("             if(N[gender_switch][jpr][apr]*ct[gender_switch][jpr][apr]*rho[gender_switch][jpr][j][apr][a]==N[g][j][a]*ct[g][j][a]*rho[g][j][jpr][a][apr]){", filename, append = TRUE, sep = "\n") #JS080416
  write("               cont_frac = 1;", filename, append = TRUE, sep = "\n") #JS080416
  write("             }", filename, append = TRUE, sep = "\n") #JS080416
  write("             else if(N[gender_switch][jpr][apr]*ct[gender_switch][jpr][apr]*rho[gender_switch][jpr][j][apr][a]!=N[g][j][a]*ct[g][j][a]*rho[g][j][jpr][a][apr] && N[g][j][a]*ct[g][j][a]*rho[g][j][jpr][a][apr] <= 0){", filename, append = TRUE, sep = "\n") #JS080416
  write("               cont_frac = 0;", filename, append = TRUE, sep = "\n") #JS080416
  write("             }", filename, append = TRUE, sep = "\n") #JS080416
  write("             else{", filename, append = TRUE, sep = "\n") #JS080416
  write("               cont_frac = (N[gender_switch][jpr][apr]*ct[gender_switch][jpr][apr]*rho[gender_switch][jpr][j][apr][a])/(N[g][j][a]*ct[g][j][a]*rho[g][j][jpr][a][apr]);", filename, append = TRUE, sep = "\n") #JS080416
  write("             }", filename, append = TRUE, sep = "\n") #JS080416
  
  write("    ", filename, append = TRUE, sep = "\n")
  write("             if(cont_frac >= 0){", filename, append = TRUE, sep = "\n")
  write("               /// fill in contact rate", filename, append = TRUE, sep = "\n")
  #write("               ctc [g][j][jpr][a][apr] = ct[gender_switch][j][apr]*pow(cont_frac, theta[g]);", filename, append = TRUE, sep = "\n") #JS080416
  write("               ctc [g][j][jpr][a][apr] = ct[g][j][a]*pow(cont_frac, theta[g]);", filename, append = TRUE, sep = "\n") #JS080416
  write("     ", filename, append = TRUE, sep = "\n")            
  write("             }", filename, append = TRUE, sep = "\n")
  write("             else{", filename, append = TRUE, sep = "\n")
  write("               ", filename, append = TRUE, sep = "\n")
  write("               if(verbose){", filename, append = TRUE, sep = "\n")
  write("                 printf(\"FOUND A NEGATIVE ENTRY IN CTC. cont_frac is: %f. REPLACING THAT BY ZERO; \", cont_frac);", filename, append = TRUE, sep = "\n")
  write("                 printf(\"N[%d][%d][%d]: %f; \", gender_switch,jpr, apr, N[gender_switch][jpr][apr]);", filename, append = TRUE, sep = "\n")
  write("                 printf(\"ct[%d][%d][%d]: %f; \", gender_switch,jpr,apr, ct[gender_switch][jpr][apr]);", filename, append = TRUE, sep = "\n")
  write("                 printf(\"rho[%d][%d][%d][%d][%d]: %f; \", gender_switch,jpr,j,apr,a, rho[gender_switch][jpr][j][apr][a]);", filename, append = TRUE, sep = "\n")
  write("                 printf(\"N[%d][%d][%d]: %f; \", g, j, a, N[g][j][a]);", filename, append = TRUE, sep = "\n")
  write("                 printf(\"ct[%d][%d][%d]: %f; \",  g, j, a, ct[g][j][a]);", filename, append = TRUE, sep = "\n")
  write("                 printf(\"rho[%d][%d][%d][%d][%d]: %f; \", g,j,jpr,a,apr,rho[g][j][jpr][a][apr]);", filename, append = TRUE, sep = "\n")
  write("               }", filename, append = TRUE, sep = "\n")
  write("            ", filename, append = TRUE, sep = "\n")
  write("               ctc [g][j][jpr][a][apr] = 0.0;", filename, append = TRUE, sep = "\n")
  write("             }", filename, append = TRUE, sep = "\n")
  write("    ", filename, append = TRUE, sep = "\n")
  write("             /// reset the contact frac ", filename, append = TRUE, sep = "\n")
  write("             cont_frac = 0.0;", filename, append = TRUE, sep = "\n")
  write("           }", filename, append = TRUE, sep = "\n")
  write("         }", filename, append = TRUE, sep = "\n")
  write("       }", filename, append = TRUE, sep = "\n")
  write("     }", filename, append = TRUE, sep = "\n")
  write("   } ", filename, append = TRUE, sep = "\n")
  
  write("     ", filename, append = TRUE, sep = "\n")
  
  write("   /// calculate lambda (force of infection)", filename, append = TRUE, sep = "\n")
  write("   double lambda[sex_length][nmb_act_classes][nmb_age_classes];", filename, append = TRUE, sep = "\n")
  write("   ", filename, append = TRUE, sep = "\n")
  write("   double inf_prob = 0.0;", filename, append = TRUE, sep = "\n")
  write("   ", filename, append = TRUE, sep = "\n")
  write("   /// loop through all combinations of traits", filename, append = TRUE, sep = "\n")
  write("   for(int g = 0; g < sex_length; g++){", filename, append = TRUE, sep = "\n")
  write("     for (int j = 0; j < nmb_act_classes; j++) {", filename, append = TRUE, sep = "\n")
  write("       for (int a = 0; a < nmb_age_classes; a++) {", filename, append = TRUE, sep = "\n")
  write("         ", filename, append = TRUE, sep = "\n")
  write("         /// define the gender switch variable (simply the opposite of the variable under consideration at the moment)", filename, append = TRUE, sep = "\n")
  write("         g == 0 ? gender_switch = 1 : gender_switch = 0;", filename, append = TRUE, sep = "\n")
  write("   ", filename, append = TRUE, sep = "\n")
  write("         /// form a sum over jpr and apr ", filename, append = TRUE, sep = "\n")
  write("         for (int jpr = 0; jpr < nmb_act_classes; jpr++) {", filename, append = TRUE, sep = "\n")
  write("           for (int apr = 0; apr < nmb_age_classes; apr++) {", filename, append = TRUE, sep = "\n")
  write("   ", filename, append = TRUE, sep = "\n")
  write("             /// ensure reasonable outcomes remain: only allow positive total populations", filename, append = TRUE, sep = "\n")
  write("             if(N [gender_switch][jpr][apr] > 0 && ctc [g][j][jpr][a][apr] > 0 && rho[g][j][jpr][a][apr] > 0 && (I_A[gender_switch][jpr][apr] + I_S[gender_switch][jpr][apr] + I_A2[gender_switch][jpr][apr] + I_S2[gender_switch][jpr][apr]) > 0){", filename, append = TRUE, sep = "\n")
  write("               inf_prob += ctc [g][j][jpr][a][apr]*rho[g][j][jpr][a][apr]*((I_A[gender_switch][jpr][apr] + I_S[gender_switch][jpr][apr] + I_A2[gender_switch][jpr][apr] + I_S2[gender_switch][jpr][apr]) / N [gender_switch][jpr][apr]);", filename, append = TRUE, sep = "\n")
  write("             }", filename, append = TRUE, sep = "\n")
  write("             else{", filename, append = TRUE, sep = "\n")
  write("               inf_prob += 0.0; ", filename, append = TRUE, sep = "\n")
  write("             }", filename, append = TRUE, sep = "\n")
  write("   ", filename, append = TRUE, sep = "\n")
  write("           }", filename, append = TRUE, sep = "\n")
  write("         }", filename, append = TRUE, sep = "\n")
  write("         lambda [g][j][a] = beta[g]*inf_prob;", filename, append = TRUE, sep = "\n") #JS010516
  write("         if(verbose){", filename, append = TRUE, sep = "\n")
  write("           /// printf(\"lambda[%d][%d][%d] is: %f; \", g,j,a, lambda [g][j][a]);", filename, append = TRUE, sep = "\n")
  write("           if(inf_prob < 0 ){", filename, append = TRUE, sep = "\n")
  write("             printf(\"FOUND A NEGATIVE ENTRY IN lambda[%d][%d][%d]. inf_prob is: %f; \", g,j,a, inf_prob);", filename, append = TRUE, sep = "\n")
  write("           }", filename, append = TRUE, sep = "\n")
  write("         }", filename, append = TRUE, sep = "\n")
  write("   ", filename, append = TRUE, sep = "\n")
  write("         /// reset to zero", filename, append = TRUE, sep = "\n")
  write("         inf_prob = 0.0;", filename, append = TRUE, sep = "\n")
  write("       }", filename, append = TRUE, sep = "\n")
  write("     }", filename, append = TRUE, sep = "\n")
  write("   } ", filename, append = TRUE, sep = "\n")
  
  write("     ", filename, append = TRUE, sep = "\n")
  write("   ///////////////////////////////////////////////////////////////////////////////", filename, append = TRUE, sep = "\n")
  write("   /// calculate the differentials (or rates of change): the results", filename, append = TRUE, sep = "\n")
  write("   ///////////////////////////////////////////////////////////////////////////////", filename, append = TRUE, sep = "\n")
  write("     ", filename, append = TRUE, sep = "\n")
  write("   double U_incoming = 0;", filename, append = TRUE, sep = "\n")
  write("   double S_incoming_age = 0;", filename, append = TRUE, sep = "\n") #JS290216
  write("   double S_incoming_vl = 0;", filename, append = TRUE, sep = "\n") #JS290216
  write("   double I_A_incoming = 0;", filename, append = TRUE, sep = "\n")
  write("   double I_S_incoming = 0;", filename, append = TRUE, sep = "\n")
  write("   double R_incoming = 0;", filename, append = TRUE, sep = "\n")
  write("   double I_A_incoming2 = 0;", filename, append = TRUE, sep = "\n")
  write("   double I_S_incoming2 = 0;", filename, append = TRUE, sep = "\n")
  write("   double Pr_stayvirginclass[sex_length][nmb_age_classes];", filename, append = TRUE, sep = "\n") #JS010316
  write("     ", filename, append = TRUE, sep = "\n")
  write("   /// outflows contains all the individuals flowing out of compartments due to age or mortality, ", filename, append = TRUE, sep = "\n")
  write("   /// and reintroduces them into the virgins compartment", filename, append = TRUE, sep = "\n")
  write("   double outflows[sex_length];", filename, append = TRUE, sep = "\n")
  write("     ", filename, append = TRUE, sep = "\n")
  write("   for(int g = 0; g < sex_length; g++){", filename, append = TRUE, sep = "\n")
  write("      outflows[g] = 0.0;", filename, append = TRUE, sep = "\n")
  write("   } ", filename, append = TRUE, sep = "\n")
  write("     ", filename, append = TRUE, sep = "\n")
  write("   /// main calcs ", filename, append = TRUE, sep = "\n")
  write("   for(int g = 0; g < sex_length; g++){", filename, append = TRUE, sep = "\n")
  write("     ", filename, append = TRUE, sep = "\n")
 
  write("        /// outflows due to aging for S, E, I_A, I_S, R, V", filename, append = TRUE, sep = "\n")
  write("        for (int jpr = 0; jpr < nmb_act_classes; jpr++) {", filename, append = TRUE, sep = "\n")
  write("           outflows[g] += aging_rate[nmb_age_classes-1]*N [g][jpr][nmb_age_classes-1];", filename, append = TRUE, sep = "\n")
  write("        }", filename, append = TRUE, sep = "\n")
  write("     ", filename, append = TRUE, sep = "\n")
  write("        /// outflows due to aging for U", filename, append = TRUE, sep = "\n")
  write("        outflows[g] += aging_rate[nmb_age_classes-1]*U[g][nmb_age_classes -1];", filename, append = TRUE, sep = "\n")
  write("        ", filename, append = TRUE, sep = "\n")

  write("     for (int a = 0; a < nmb_age_classes; a++) {", filename, append = TRUE, sep = "\n") #JS010316
  write("        if(a > 0){ ", filename, append = TRUE, sep = "\n") #JS010316
  write("           if(virg_loss[g][a-1]>0){", filename, append = TRUE, sep = "\n") #JS010316
  write("              Pr_stayvirginclass[g][a] = virg_loss[g][a] / virg_loss[g][a-1];", filename, append = TRUE, sep = "\n") #JS010316
  write("           }", filename, append = TRUE, sep = "\n") #JS010316
  write("           else {", filename, append = TRUE, sep = "\n") #JS010316
  write("              Pr_stayvirginclass[g][a] = 1;", filename, append = TRUE, sep = "\n") #JS010316
  write("           }", filename, append = TRUE, sep = "\n") #JS010316
  write("         }", filename, append = TRUE, sep = "\n") #JS010316
  write("         else {", filename, append = TRUE, sep = "\n") #if a==0
  write("           if(virg_loss[g][a]>0){", filename, append = TRUE, sep = "\n") #JS010316
  write("              Pr_stayvirginclass[g][a] = virg_loss[g][a];", filename, append = TRUE, sep = "\n") #JS010316
  write("           }", filename, append = TRUE, sep = "\n") #JS010316
  write("           else {", filename, append = TRUE, sep = "\n") #JS010316
  write("              Pr_stayvirginclass[g][a] = 0;", filename, append = TRUE, sep = "\n") #JS010316
  write("           }", filename, append = TRUE, sep = "\n") #JS010316
  write("         }", filename, append = TRUE, sep = "\n") #JS010316
  write("     }", filename, append = TRUE, sep = "\n") #JS010316
  
  write("    ", filename, append = TRUE, sep = "\n")
  write("     for (int j = 0; j < nmb_act_classes; j++) {", filename, append = TRUE, sep = "\n")
  write("       for (int a = 0; a < nmb_age_classes; a++) {", filename, append = TRUE, sep = "\n")
  write("    ", filename, append = TRUE, sep = "\n")
  
  write("         /// only assign variable once for U (there always must be one activity class, j =0, so the assignment will always be performed)", filename, append = TRUE, sep = "\n")
  write("         if (j == 0) {", filename, append = TRUE, sep = "\n")
  write("           ", filename, append = TRUE, sep = "\n")
  write("           if(a >= 1){", filename, append = TRUE, sep = "\n")
  write("             U_incoming = aging_rate[a-1] * Pr_stayvirginclass[g][a] * U [g][a-1];", filename, append = TRUE, sep = "\n") #JS290216
  write("           }", filename, append = TRUE, sep = "\n")
  write("           else{", filename, append = TRUE, sep = "\n")
  write("             U_incoming = Pr_stayvirginclass[g][a] * outflows[g];", filename, append = TRUE, sep = "\n") #JS030316
  write("           }", filename, append = TRUE, sep = "\n")
  write("     ", filename, append = TRUE, sep = "\n")
  write("           dU [g][a] = - aging_rate[a] * U [g][a] + U_incoming;", filename, append = TRUE, sep = "\n") #290216
  write("         }", filename, append = TRUE, sep = "\n")
  write("    ", filename, append = TRUE, sep = "\n")
  
  write("         /// number of incoming individuals is zero if a is = 0, and corresponds to the value of the one-unit-younger age group", filename, append = TRUE, sep = "\n")
  write("         if (a >= 1) {", filename, append = TRUE, sep = "\n")
  write("           /// the incoming stem from the cohort which is one \"a\"-unit lower than the one we are looping through", filename, append = TRUE, sep = "\n")
  write("           S_incoming_age = aging_rate[a-1] * S [g][j][a - 1];", filename, append = TRUE, sep = "\n") #290216
  write("           S_incoming_vl = aging_rate[a-1] * (1-Pr_stayvirginclass[g][a]) * U [g][a - 1];", filename, append = TRUE, sep = "\n") #290216
  write("           I_A_incoming = aging_rate[a-1] * I_A [g][j][a - 1];", filename, append = TRUE, sep = "\n")
  write("           I_S_incoming = aging_rate[a-1] * I_S [g][j][a - 1];", filename, append = TRUE, sep = "\n")
  write("           R_incoming = aging_rate[a-1] * R [g][j][a - 1];", filename, append = TRUE, sep = "\n")
  write("           I_A_incoming2 = aging_rate[a-1] * I_A2 [g][j][a - 1];", filename, append = TRUE, sep = "\n")
  write("           I_S_incoming2 = aging_rate[a-1] * I_S2 [g][j][a - 1];", filename, append = TRUE, sep = "\n")
  write("         } ", filename, append = TRUE, sep = "\n")
  write("         else {", filename, append = TRUE, sep = "\n")
  write("           /// in the smallest cohort, there are no incoming except for the ones that are born. ", filename, append = TRUE, sep = "\n")
  write("           S_incoming_age = 0.0;", filename, append = TRUE, sep = "\n") #JS290216
  write("           S_incoming_vl = (1-Pr_stayvirginclass[g][a]) * outflows[g];", filename, append = TRUE, sep = "\n") #JS030316
  write("           I_A_incoming = 0.0;", filename, append = TRUE, sep = "\n")
  write("           I_S_incoming = 0.0;", filename, append = TRUE, sep = "\n")
  write("           R_incoming = 0.0;", filename, append = TRUE, sep = "\n")
  write("           I_A_incoming2 = 0.0;", filename, append = TRUE, sep = "\n")
  write("           I_S_incoming2 = 0.0;", filename, append = TRUE, sep = "\n")
  write("         }", filename, append = TRUE, sep = "\n")
  write("     ", filename, append = TRUE, sep = "\n")
  
  write("         /// calculate the differentials ", filename, append = TRUE, sep = "\n")
  write("         double tottreatments = screen[g][a] * (I_A[g][j][a]+I_A2[g][j][a]) + treat * (I_S[g][j][a]+I_S2[g][j][a]);", filename, append = TRUE, sep = "\n") 
  write("         dS [g][j][a] = -aging_rate[a]*S[g][j][a] + fold[g][j] * S_incoming_vl + S_incoming_age - lambda[g][j][a]*S[g][j][a] + (omega[0] * screen[g][a]*I_A[g][j][a] + omega[1] * treat * I_S[g][j][a]) ;", filename, append = TRUE, sep = "\n") #JS290216
  write("         dI_A [g][j][a] = -aging_rate[a]*I_A[g][j][a] + I_A_incoming + (1-fsymp[g])* lambda[g][j][a]*S[g][j][a] - (gamma[g] + omega[0] * screen[g][a])*I_A[g][j][a] ;", filename, append = TRUE, sep = "\n")
  write("         dI_S [g][j][a] = -aging_rate[a]*I_S[g][j][a] + I_S_incoming + fsymp[g] * lambda[g][j][a]*S[g][j][a] - (omega[1] * treat)*I_S[g][j][a];", filename, append = TRUE, sep = "\n")
  write("         dR [g][j][a] = -aging_rate[a]*R[g][j][a] + R_incoming + gamma[g]*I_A[g][j][a]+ gamma[g]*I_A2[g][j][a] - ((1-kappa[g])* lambda[g][j][a])*R[g][j][a] + (omega[0] * screen[g][a]*I_A2[g][j][a] + omega[1] * treat * I_S2[g][j][a]) ;", filename, append = TRUE, sep = "\n")
  write("         dI_A2 [g][j][a] = -aging_rate[a]*I_A2[g][j][a] + I_A_incoming2 + (1-fsymp[g])* (1-kappa[g]) * lambda[g][j][a]*R[g][j][a] - (gamma[g] + omega[0] * screen[g][a])*I_A2[g][j][a] ;", filename, append = TRUE, sep = "\n")
  write("         dI_S2 [g][j][a] = -aging_rate[a]*I_S2[g][j][a] + I_S_incoming2 + fsymp[g] * (1-kappa[g]) * lambda[g][j][a]*R[g][j][a] - (omega[1] * treat)*I_S2[g][j][a];", filename, append = TRUE, sep = "\n")
  write("         dD [g][j][a] = tottreatments;", filename, append = TRUE, sep = "\n") # treatment efficacy not used here: also unsuccesful treatments get in the testing records
  write("         dInc [g][j][a] = lambda[g][j][a]*S[g][j][a] + (1-kappa[g]) * lambda[g][j][a]*R[g][j][a];", filename, append = TRUE, sep = "\n") # treatment efficacy not used here: also unsuccesful treatments get in the testing records
  
  write("     ", filename, append = TRUE, sep = "\n")
  write("         /// switching between classes ", filename, append = TRUE, sep = "\n")
  write("         for (int jpr = 0; jpr < nmb_act_classes; jpr++) {", filename, append = TRUE, sep = "\n")
  write("     ", filename, append = TRUE, sep = "\n")
  write("           dS[g][j][a] -= sw [g][j][jpr][a] * S [g][j][a];", filename, append = TRUE, sep = "\n")
  write("           dS[g][j][a] += sw [g][jpr][j][a] * S [g][jpr][a];", filename, append = TRUE, sep = "\n")
  write("     ", filename, append = TRUE, sep = "\n")
  write("           dI_A[g][j][a] -= sw [g][j][jpr][a] * I_A [g][j][a];", filename, append = TRUE, sep = "\n")
  write("           dI_A[g][j][a] += sw [g][jpr][j][a] * I_A [g][jpr][a];", filename, append = TRUE, sep = "\n")
  write("     ", filename, append = TRUE, sep = "\n")
  write("           dI_S[g][j][a] -= sw [g][j][jpr][a] * I_S [g][j][a];", filename, append = TRUE, sep = "\n")
  write("           dI_S[g][j][a] += sw [g][jpr][j][a] * I_S [g][jpr][a];", filename, append = TRUE, sep = "\n")
  write("     ", filename, append = TRUE, sep = "\n")
  write("           dR[g][j][a] -= sw [g][j][jpr][a] * R [g][j][a];", filename, append = TRUE, sep = "\n")
  write("           dR[g][j][a] += sw [g][jpr][j][a] * R [g][jpr][a];", filename, append = TRUE, sep = "\n")
  write("     ", filename, append = TRUE, sep = "\n")
  write("           dI_A2[g][j][a] -= sw [g][j][jpr][a] * I_A2 [g][j][a];", filename, append = TRUE, sep = "\n")
  write("           dI_A2[g][j][a] += sw [g][jpr][j][a] * I_A2 [g][jpr][a];", filename, append = TRUE, sep = "\n")
  write("     ", filename, append = TRUE, sep = "\n")
  write("           dI_S2[g][j][a] -= sw [g][j][jpr][a] * I_S2 [g][j][a];", filename, append = TRUE, sep = "\n")
  write("           dI_S2[g][j][a] += sw [g][jpr][j][a] * I_S2 [g][jpr][a];", filename, append = TRUE, sep = "\n")
  write("           ", filename, append = TRUE, sep = "\n")
  write("         }", filename, append = TRUE, sep = "\n")
  write("     ", filename, append = TRUE, sep = "\n")      
  write("         /// reset everything", filename, append = TRUE, sep = "\n")
  write("         U_incoming = 0.0;", filename, append = TRUE, sep = "\n")
  write("         S_incoming_age = 0.0;", filename, append = TRUE, sep = "\n") #JS290216
  write("         S_incoming_vl = 0.0;", filename, append = TRUE, sep = "\n") #JS290216
  write("         I_A_incoming = 0.0;", filename, append = TRUE, sep = "\n")
  write("         I_S_incoming = 0.0;", filename, append = TRUE, sep = "\n")
  write("         R_incoming = 0.0;", filename, append = TRUE, sep = "\n")
  write("         I_A_incoming2 = 0.0;", filename, append = TRUE, sep = "\n")
  write("         I_S_incoming2 = 0.0;", filename, append = TRUE, sep = "\n")
  write("         ", filename, append = TRUE, sep = "\n")
  write("       }", filename, append = TRUE, sep = "\n")
  write("     }", filename, append = TRUE, sep = "\n")
  write("   } ", filename, append = TRUE, sep = "\n")
  
  
  write("   ///////////////////////////////////////////////////////////////////////////////", filename, append = TRUE, sep = "\n")
  write("     /// return the results", filename, append = TRUE, sep = "\n")
  write("   ///////////////////////////////////////////////////////////////////////////////", filename, append = TRUE, sep = "\n")
  
  write("     ", filename, append = TRUE, sep = "\n")
  write("   int res_counter = 0; ", filename, append = TRUE, sep = "\n")
  
  write("     ", filename, append = TRUE, sep = "\n")
  write("   for(int g = 0; g < sex_length; g++){", filename, append = TRUE, sep = "\n")
  write("     for (int a = 0; a < nmb_age_classes; a++) {", filename, append = TRUE, sep = "\n")
  write("       res [res_counter] = dU [g][a];", filename, append = TRUE, sep = "\n")
  write("       res_counter++;", filename, append = TRUE, sep = "\n")
  write("     }", filename, append = TRUE, sep = "\n")
  write("   }", filename, append = TRUE, sep = "\n")
  write("     ", filename, append = TRUE, sep = "\n")
  write("   for(int g = 0; g < sex_length; g++){", filename, append = TRUE, sep = "\n")
  write("     for (int j = 0; j < nmb_act_classes; j++) {", filename, append = TRUE, sep = "\n")
  write("       for (int a = 0; a < nmb_age_classes; a++) {", filename, append = TRUE, sep = "\n")
  write("         res [res_counter] = dS [g][j][a];", filename, append = TRUE, sep = "\n")
  write("         res [res_counter + 1*tensor_entry_numbers] = dI_A [g][j][a];", filename, append = TRUE, sep = "\n")
  write("         res [res_counter + 2*tensor_entry_numbers] = dI_S [g][j][a];", filename, append = TRUE, sep = "\n")
  write("         res [res_counter + 3*tensor_entry_numbers] = dR [g][j][a];", filename, append = TRUE, sep = "\n")
  write("         res [res_counter + 4*tensor_entry_numbers] = dI_A2 [g][j][a];", filename, append = TRUE, sep = "\n")
  write("         res [res_counter + 5*tensor_entry_numbers] = dI_S2 [g][j][a];", filename, append = TRUE, sep = "\n")
  write("         res [res_counter + 6*tensor_entry_numbers] = dD [g][j][a];", filename, append = TRUE, sep = "\n")
  write("         res [res_counter + 7*tensor_entry_numbers] = dInc[g][j][a];", filename, append = TRUE, sep = "\n")
  write("   ", filename, append = TRUE, sep = "\n")
  write("         res_counter++;", filename, append = TRUE, sep = "\n")
  write("       }", filename, append = TRUE, sep = "\n")
  write("     }", filename, append = TRUE, sep = "\n")
  write("   }", filename, append = TRUE, sep = "\n")
  
#   write("     ", filename, append = TRUE, sep = "\n")
#   write("   if (verbose) {", filename, append = TRUE, sep = "\n")
#   write("     printf(\"######################## Exiting the program ################################################\");", filename, append = TRUE, sep = "\n")
#   write("   }", filename, append = TRUE, sep = "\n")
  
  #########################################################################################################################
  ### end of the function
  write(c("   }","}"), filename, sep = "\n", append = TRUE)
  
}

########## modeling function 
### this is a function that calculates the differentials, given some state variable values (y) and parameters (parms)
SIR.cpp <- function(t, y.cpp, parms.cpp){
  
  nsex <- 2
  nAC <- (length(y.cpp)-nsex)/(4*nJ*nsex+nsex)
  
  ### convert parameters into a string
  parms.string <- paste(parms.cpp, sep = " ")
  init.string <- paste(y.cpp, sep = " ")
  
  ### combine strings into one arguments string vector
  arg.string <- paste(c(t,parms.string, init.string)) #JS150316
  
  ### declare a results vector of the right length
  res <- numeric(length(init.string))
  
  ### get out the differentials
  out <- .C("STImod",as.double(arg.string),as.double(res))[[2]]
  
  ##### out goes:
  return(list(res = out))
}

SIR.cpp <- cmpfun(SIR.cpp)
#res <- SIR.cpp(t = 0, init.cpp, parameters.cpp)