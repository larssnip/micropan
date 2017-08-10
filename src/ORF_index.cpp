

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::DataFrame ORF_index(SEXP Tags, SEXP Sequence) {
  std::vector<std::string> tag = as<std::vector<std::string> >(Tags);
  std::vector<std::string> seq  = as<std::vector<std::string> >(Sequence);
  std::vector<int> start1(10000,0);
  std::vector<int> start2(10000,0);
  std::vector<int> start3(10000,0);
  std::vector<int> start1rc(10000,0);
  std::vector<int> start2rc(10000,0);
  std::vector<int> start3rc(10000,0);
  int s1=-1,   s2=-1,   s3=-1, e1=-1,   e2=-1,   e3=-1;
  bool st1=true, st2=true, st3=true;
  int n = seq.size();
  
  std::vector<std::vector<std::string> > out; // Extracted sequences
  std::vector<std::string> outGenomeSequence;
  std::vector<int> outStrand;
  std::vector<int> outLeft;
  std::vector<int> outRight;
  std::vector<int> outPartial;
  
  // Loop over Genome Sequences
  for(int str=0; str<n; ++str){
    std::vector<int> outoutLeft;
    std::vector<int> outoutRight;
    std::vector<int> outoutLeftRC;
    std::vector<int> outoutRightRC;
    std::vector<int> outoutPartial;
    std::vector<int> outoutPartialRC;
    
    std::string cur = seq[str];
    int lStr = cur.length();
    int ret = 0;
    s1=-1; s2=-1; s3=-1; e1=-1; e2=-1; e3=-1;
    st1=true; st2=true; st3=true;

    
    std::string C1 = cur.substr(0,1);
    std::string C2 = cur.substr(0,1);
    std::string C3 = cur.substr(1,1);
    
    // Loop over positions in Genome Sequence
    for(int pos=2; pos <lStr; ++pos){
      
      // Iterate position
      C1 = C2;
      C2 = C3;
      C3 = cur[pos];
      
      // Minimal testing for start/end forward/reverse complement (start = 1/-1, stop = 2/-2)
      ret = 0;
      if(C1 == "T"){
        if(C3 == "G"){
          if(C2 == "A"){
            ret = 2;
          } else {
            if(C2 == "T"){
              ret = 1;
            }
          }
        } else {
          if(C3 == "A"){
            if(C2 == "A"){
              ret = 2;
            } else {
              if(C2 == "G"){
                ret = 2;
              } else {
                if(C2 == "C"){
                  ret = -2;
                } else {
                  if(C2 == "T"){
                    ret = -2;
                  }
                }
              }
            }
          }
        }
      } else {
        if(C1 == "C"){
          if(C2 == "A"){
            if(C3 != "G"){
              ret = -1;
            }
          } else {
            if(C2 == "T"){
              if(C3 == "A"){
                ret = -2;
              }
            }
          }
        } else {
          if(C2 == "T"){
            if(C3 == "G"){
              if(C1 == "A"){
                ret = 1;
              } else {
                if(C1 == "G"){
                  ret = 1;
                }
              }
            }
          }
        }
      } // end of test sequence
      
      // Take consequences of the ret value (-2,-1,0,1,2)
      if(ret != 0){
        if(ret > 0){ 
          if(ret == 1){ // Found start codon on forward strand
            if(pos % 3 == 2){ // s1, start1
              ++s1;
              start1[s1] = pos;
            } else {
              if(pos % 3 == 0){ // s2, start2
                ++s2;
                start2[s2] = pos;
              } else { // s3, start3
                ++s3;
                start3[s3] = pos;
              }
            }
          } else { // ret == 2, Found stop codon on forward strand
            if(pos % 3 == 2){ // s1, start1, first reading frame
              if(s1 == -1 && st1){ // Partial match, no start codon forward
                outoutLeft.push_back( 1 );
                outoutRight.push_back( pos+1 );
                outoutPartial.push_back( 1 );
                st1 = false;
              } else { // Full match forward
                for(int ss1=0; ss1<=s1; ++ss1){
                  outoutLeft.push_back( start1[ss1]-1 );
                  outoutRight.push_back( pos+1 );
                  outoutPartial.push_back( 0 );
                }
                s1  = -1;
                st1 = false;
              }
            } else {
              if(pos % 3 == 0){ // s2, start2
                if(s2 == -1 && st2){ // Partial match, no start codon forward
                  outoutLeft.push_back( 2 );
                  outoutRight.push_back( pos+1 );
                  outoutPartial.push_back( 1 );
                  st2 = false;
                } else { // Full match forward
                  for(int ss2=0; ss2<=s2; ++ss2){
                    outoutLeft.push_back( start2[ss2]-1 );
                    outoutRight.push_back( pos+1 );
                    outoutPartial.push_back( 0 );
                  }
                  s2  = -1;
                  st2 = false;
                }
              } else { // s3, start3
                if(s3 == -1 && st3){ // Partial match, no start codon forward
                  outoutLeft.push_back( 3 );
                  outoutRight.push_back( pos+1 );
                  outoutPartial.push_back( 1 );
                  st3 = false;
                } else { // Full match forward
                  for(int ss3=0; ss3<=s3; ++ss3){
                    outoutLeft.push_back( start3[ss3]-1 );
                    outoutRight.push_back( pos+1 );
                    outoutPartial.push_back( 0 );
                  }
                  s3  = -1;
                  st3 = false;
                }
              }
            }
          }
        } else { // Reverse complement  --<--[---------<---<-------<----------[----<---<--<----
          if(ret == -1){ // Combine with endpoint
            if(pos % 3 == 2){ // e1
              if(e1 < 0){ // No stop codon
                outoutLeftRC.push_back( 1 );
                outoutRightRC.push_back( pos+1 );
                outoutPartialRC.push_back( 1 );
              } else { // Full match
                outoutLeftRC.push_back( e1-1 );
                outoutRightRC.push_back( pos+1 );
                outoutPartialRC.push_back( 0 );
              }
            } else {
              if(pos % 3 == 0){ // e2
                if(e2 < 0){ // No stop codon
                  outoutLeftRC.push_back( 2 );
                  outoutRightRC.push_back( pos+1 );
                  outoutPartialRC.push_back( 1 );
                } else { // Full match
                  outoutLeftRC.push_back( e2-1 );
                  outoutRightRC.push_back( pos+1 );
                  outoutPartialRC.push_back( 0 );
                }
              } else { // s3, start3
                if(e3 < 0){ // No stop codon
                  outoutLeftRC.push_back( 3 );
                  outoutRightRC.push_back( pos+1 );
                  outoutPartialRC.push_back( 1 );
                } else { // Full match
                  outoutLeftRC.push_back( e3-1 );
                  outoutRightRC.push_back( pos+1 );
                  outoutPartialRC.push_back( 0 );
                }
              }
            }
          } else { // ret == -2, store endpoint
            if(pos % 3 == 2){ // s1, start1
              e1 = pos;
            } else {
              if(pos % 3 == 0){ // s2, start2
                e2 = pos;
              } else { // s3, start3
                e3 = pos;
              }
            }
          }
        }
        
        // Reset ret
        ret = 0;
      } // ret != 0

      if(pos >=lStr-3){ // Last letter of last codon
        if(ret <= 0){ // Last letter does not make start/stop codon (forward)
          if(pos % 3 == 2 && s1 >= 0){
            for(int ss1=0; ss1<=s1; ++ss1){
              outoutLeft.push_back( start1[ss1]-1 );
              outoutRight.push_back( lStr-2 );
              outoutPartial.push_back( -1 );
            }
          }
          if(pos % 3 == 0 && s2 >= 0){
            for(int ss2=0; ss2<=s2; ++ss2){
              outoutLeft.push_back( start2[ss2]-1 );
              outoutRight.push_back( lStr-1 );
              outoutPartial.push_back( -1 );
            }
          }
          if(pos % 3 == 1 && s3 >= 0){
            for(int ss3=0; ss3<=s3; ++ss3){
              outoutLeft.push_back( start3[ss3]-1 );
              outoutRight.push_back( lStr );
              outoutPartial.push_back( -1 );
            }
          }
        }
        if(ret >= 0){ // Last letter does not make start/stop codon (reverse complement)
          if(pos % 3 == 2 && e1 >= 0){
            outoutLeftRC.push_back( e1-1 );
            outoutRightRC.push_back( lStr-2 );
            outoutPartialRC.push_back( -1 );
          }
          if(pos % 3 == 0 && e2 >= 0){
            outoutLeftRC.push_back( e2-1 );
            outoutRightRC.push_back( lStr-1 );
            outoutPartialRC.push_back( -1 );
          }
          if(pos % 3 == 1 && e3 >= 0){
            outoutLeftRC.push_back( e3-1 );
            outoutRightRC.push_back( lStr );
            outoutPartialRC.push_back( -1 );
          }
        }
      } // Last three
    } // end loop over positions in Genome Sequence
    
    
    // Store results in main vectors
    outLeft.insert(  outLeft.end(),  outoutLeft.begin(),  outoutLeft.end()  );
    outRight.insert( outRight.end(), outoutRight.begin(), outoutRight.end() );
    outLeft.insert(  outLeft.end(),  outoutLeftRC.begin(),  outoutLeftRC.end()  );
    outRight.insert( outRight.end(), outoutRightRC.begin(), outoutRightRC.end() );
    std::vector<int> outoutStrandP(outoutLeft.size(),1);
    std::vector<int> outoutStrandM(outoutLeftRC.size(),-1);
    outStrand.insert( outStrand.end(), outoutStrandP.begin(), outoutStrandP.end() );
    outStrand.insert( outStrand.end(), outoutStrandM.begin(), outoutStrandM.end() );
    std::vector<std::string> outoutGenomeSequence(outoutLeft.size()+outoutLeftRC.size(),tag[str]);
    outGenomeSequence.insert( outGenomeSequence.end(), outoutGenomeSequence.begin(), outoutGenomeSequence.end() );
    outPartial.insert(  outPartial.end(),  outoutPartial.begin(),  outoutPartial.end()  );
    outPartial.insert(  outPartial.end(),  outoutPartialRC.begin(),  outoutPartialRC.end()  );
    
  } // end loop over Genome Sequences
  
  
  // Return av data
  Rcpp::DataFrame outDF;
  outDF =
  Rcpp::DataFrame::create(
    Rcpp::Named("GenomeSequence") = outGenomeSequence,
    Rcpp::Named("Strand")         = outStrand,
    Rcpp::Named("Left")           = outLeft,
    Rcpp::Named("Right")          = outRight,
    Rcpp::Named("Partial")        = outPartial);
    
  return outDF;
}
