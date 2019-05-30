#include <math.h>
#include <Rmath.h>
#include <stdlib.h>
#include <R.h>
#include <stdbool.h>
#include <stdio.h>

// Quickly run through RL trials
void updateValuesC(int *nTrials, 
                   double *value1, double *value2,
                   double *VV1, double *VV2, 
                   int *choice, int *outcome, 
                   double *PE, 
                   double *eta1, 
                   double *eta2) {
  
  // declare some variables
  double updateValue = 0.;
  double this_eta1 = 0;
  double this_eta2 = 0;
  int c_ = 0;
  int o_ = 0;
  double dv = 0.;
  double v1 = *value1;
  double v2 = *value2;
  
  for(unsigned int i = 0; i < *nTrials; i++) {
    c_ = choice[i]; // c_ = choice
    o_ = outcome[i]; // o_ = outcome
    this_eta1 = eta1[i];  // trialwise learning rate for positive PE
    this_eta2 = eta2[i];  // trialwise learning rate for negative PE
    
    // Calculate PE (="dv" for delta value)
    if(c_ == 1) {
      dv = o_ - v1;
    } else {
      dv = o_ - v2;
    }
    
    // output
    VV1[i] = v1;
    VV2[i] = v2;
    PE[i] = dv;
    
    // get learning rate
    if(dv > 0) {
      updateValue = this_eta1;
    } else {
      updateValue = this_eta2;
    }
    
    // Update values
    if(c_ == 1) {
      v1 = v1 + updateValue*dv;
    } else {
      v2 = v2 + updateValue*dv;
    }
  }
}


// Quickly run through RL trials
void doubleUpdateValuesC(int *nTrials, 
                         double *value1, double *value2,
                         double *VV1, double *VV2, 
                         int *choice, int *outcome1, int *outcome2,
                         double *PE1, double *PE2,
                         double *eta1, double *eta2) {
  
  // declare some variables
  double updateValue1 = 0.;
  double updateValue2 = 0.;
  int o1_ = 0;
  int o2_ = 0;
  double dv1 = 0.;
  double dv2 = 0.;
  double v1 = *value1;
  double v2 = *value2;
  
  for(unsigned int i = 0; i < *nTrials; i++) {
    o1_ = outcome1[i]; // o1_ = outcome1
    o2_ = outcome2[i]; // o2_ = outcome2
    
    // Calculate PE (dv here)
    dv1 = o1_ - v1;
    dv2 = o2_ - v2;

    // output
    VV1[i] = v1;
    VV2[i] = v2;
    PE1[i] = dv1;
    PE2[i] = dv2;
    
    // get learning rate
    if(dv1 > 0) {
      updateValue1 = *eta1;
    } else {
      updateValue1 = *eta2;
    }
    
    if(dv2 > 0) {
      updateValue2 = *eta1;
    } else {
      updateValue2 = *eta2;
    }
    
    // Update values
    v1 = v1 + updateValue1*dv1;
    v2 = v2 + updateValue2*dv2;
  }
}

/*
void tester(double *VV, int *nrow, int *ncol) {
  
  //printf("Value 0: %.3f", VV[0]);
  
  int r_loc = *nrow;
  int r_col = *ncol;
  for(int row=0; row < r_loc; row++) {
    printf("Row %d", row);
    for(int col=0; col < r_col; col++) {
      printf(", col %d",col);
      if(VV[col*r_loc + row] != VV[col*r_loc + row]) {
        printf("Value is NA!!!!!!!!!!!! WOOOT");
      } else {
        printf(", value: %.3f\n", VV[col*r_loc + row]);
      }
    }
  }
}
*/

/// CMC = updateValues in C for multiple choice
void doubleUpdateValuesCMC(int *nTrials, int *nChoices,
                           double *values,
                           double *VV, // of length ntrials*n_pairs
                           double *PE,
                           double *outcome,
                           double *eta1, 
                           double *eta2) {
  // nTrials: number of trials
  // nChoices: total number of choice options. For example, 3 pairs of 2 stimuli = 6 choices
  // values: initial value of each choice option (length nChoices)
  // VV: output for trial-by-trial values of each choice option. In R, a matrix of size (nTrials, nChoices); here in C, an array of length (nTrials*nChoices)
  // PE: output for trial-by-trial prediction errors for each choie option. Identical size as VV
  // outcome: (nTrials, nChoices): trial by trial outcomes per choice option. NA is no output.
  // eta1, eta2: floats with the learning rates (positive and negative, respectively)
  
  // declare some variables
  double updateValue;
  double this_eta1 = 0;
  double this_eta2 = 0;
  int mat_idx = 0;
  int nt = *nTrials;
  int nc = *nChoices;
  static double dv[64] = {0};   // WARNING: MAXIMUM NUMBER OF CHOICES IS FIXED HERE TO 64. Anything more will crash.

  // Loop over trials i
  for(unsigned int i = 0; i < nt; i++) {
//    printf("Trial N: %d\n", i);
    this_eta1 = eta1[i];  // trialwise learning rate for positive PE
    this_eta2 = eta2[i];  // trialwise learning rate for negative PE
    
    // Loop over choice alternatives
    for(unsigned int ch = 0; ch < nc; ch++) {
//      printf("Choice option: %d\n", ch);
      mat_idx = nt*ch+i;  // Where in the VV-matrix and outcome-matrix are we?
      VV[mat_idx] = values[ch]; // Offload current values
      
      // If the outcome for this choice was not NA, update (note that o_ == o_ returns FALSE if o_ is NA)
      if(outcome[mat_idx] == outcome[mat_idx]) {
//        printf("Outcome is: %.3f, this is NOT NA!\n", outcome[mat_idx]);
        dv[ch] = outcome[mat_idx] - values[ch];  // prediction error dv = outcome - value
        PE[mat_idx] = dv[ch];  // offload PE
        
        // Do we update with eta1 or eta2?
        updateValue = dv[ch] > 0 ? this_eta1 : this_eta2;  // Ternary expression  (if-else in a one-liner)
  
        values[ch] = values[ch] + updateValue*dv[ch];  // Update value
      }
    }
  }
}