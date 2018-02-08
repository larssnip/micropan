

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
  std::vector<std::string> outSeqid;
  std::vector<int> outStrand;
  std::vector<int> outStart;
  std::vector<int> outEnd;
  std::vector<int> outTruncated;
  
  // Loop over Genome Sequences
  for(int str=0; str<n; ++str){
    std::vector<int> outoutStart;
    std::vector<int> outoutEnd;
    std::vector<int> outoutStartRC;
    std::vector<int> outoutEndRC;
    std::vector<int> outoutTruncated;
    std::vector<int> outoutTruncatedRC;
    
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
      
      // Minimal testing for start/end forward/reverse complement (start = 1/-1, End = 2/-2)
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
              if(s1 == -1 && st1){ // Truncated match, no start codon forward
                outoutStart.push_back( 1 );
                outoutEnd.push_back( pos+1 );
                outoutTruncated.push_back( 1 );
                st1 = false;
              } else { // Full match forward
                for(int ss1=0; ss1<=s1; ++ss1){
                  outoutStart.push_back( start1[ss1]-1 );
                  outoutEnd.push_back( pos+1 );
                  outoutTruncated.push_back( 0 );
                }
                s1  = -1;
                st1 = false;
              }
            } else {
              if(pos % 3 == 0){ // s2, start2
                if(s2 == -1 && st2){ // Truncated match, no start codon forward
                  outoutStart.push_back( 2 );
                  outoutEnd.push_back( pos+1 );
                  outoutTruncated.push_back( 1 );
                  st2 = false;
                } else { // Full match forward
                  for(int ss2=0; ss2<=s2; ++ss2){
                    outoutStart.push_back( start2[ss2]-1 );
                    outoutEnd.push_back( pos+1 );
                    outoutTruncated.push_back( 0 );
                  }
                  s2  = -1;
                  st2 = false;
                }
              } else { // s3, start3
                if(s3 == -1 && st3){ // Truncated match, no start codon forward
                  outoutStart.push_back( 3 );
                  outoutEnd.push_back( pos+1 );
                  outoutTruncated.push_back( 1 );
                  st3 = false;
                } else { // Full match forward
                  for(int ss3=0; ss3<=s3; ++ss3){
                    outoutStart.push_back( start3[ss3]-1 );
                    outoutEnd.push_back( pos+1 );
                    outoutTruncated.push_back( 0 );
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
                outoutStartRC.push_back( 1 );
                outoutEndRC.push_back( pos+1 );
                outoutTruncatedRC.push_back( 1 );
              } else { // Full match
                outoutStartRC.push_back( e1-1 );
                outoutEndRC.push_back( pos+1 );
                outoutTruncatedRC.push_back( 0 );
              }
            } else {
              if(pos % 3 == 0){ // e2
                if(e2 < 0){ // No stop codon
                  outoutStartRC.push_back( 2 );
                  outoutEndRC.push_back( pos+1 );
                  outoutTruncatedRC.push_back( 1 );
                } else { // Full match
                  outoutStartRC.push_back( e2-1 );
                  outoutEndRC.push_back( pos+1 );
                  outoutTruncatedRC.push_back( 0 );
                }
              } else { // s3, start3
                if(e3 < 0){ // No stop codon
                  outoutStartRC.push_back( 3 );
                  outoutEndRC.push_back( pos+1 );
                  outoutTruncatedRC.push_back( 1 );
                } else { // Full match
                  outoutStartRC.push_back( e3-1 );
                  outoutEndRC.push_back( pos+1 );
                  outoutTruncatedRC.push_back( 0 );
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
              outoutStart.push_back( start1[ss1]-1 );
              outoutEnd.push_back( lStr-2 );
              outoutTruncated.push_back( -1 );
            }
          }
          if(pos % 3 == 0 && s2 >= 0){
            for(int ss2=0; ss2<=s2; ++ss2){
              outoutStart.push_back( start2[ss2]-1 );
              outoutEnd.push_back( lStr-1 );
              outoutTruncated.push_back( -1 );
            }
          }
          if(pos % 3 == 1 && s3 >= 0){
            for(int ss3=0; ss3<=s3; ++ss3){
              outoutStart.push_back( start3[ss3]-1 );
              outoutEnd.push_back( lStr );
              outoutTruncated.push_back( -1 );
            }
          }
        }
        if(ret >= 0){ // Last letter does not make start/stop codon (reverse complement)
          if(pos % 3 == 2 && e1 >= 0){
            outoutStartRC.push_back( e1-1 );
            outoutEndRC.push_back( lStr-2 );
            outoutTruncatedRC.push_back( -1 );
          }
          if(pos % 3 == 0 && e2 >= 0){
            outoutStartRC.push_back( e2-1 );
            outoutEndRC.push_back( lStr-1 );
            outoutTruncatedRC.push_back( -1 );
          }
          if(pos % 3 == 1 && e3 >= 0){
            outoutStartRC.push_back( e3-1 );
            outoutEndRC.push_back( lStr );
            outoutTruncatedRC.push_back( -1 );
          }
        }
      } // Last three
    } // end loop over positions in Genome Sequence
    
    
    // Store results in main vectors
    outStart.insert(  outStart.end(),  outoutStart.begin(),  outoutStart.end()  );
    outEnd.insert( outEnd.end(), outoutEnd.begin(), outoutEnd.end() );
    outStart.insert(  outStart.end(),  outoutStartRC.begin(),  outoutStartRC.end()  );
    outEnd.insert( outEnd.end(), outoutEndRC.begin(), outoutEndRC.end() );
    std::vector<int> outoutStrandP(outoutStart.size(),1);
    std::vector<int> outoutStrandM(outoutStartRC.size(),-1);
    outStrand.insert( outStrand.end(), outoutStrandP.begin(), outoutStrandP.end() );
    outStrand.insert( outStrand.end(), outoutStrandM.begin(), outoutStrandM.end() );
    std::vector<std::string> outoutSeqid(outoutStart.size()+outoutStartRC.size(),tag[str]);
    outSeqid.insert( outSeqid.end(), outoutSeqid.begin(), outoutSeqid.end() );
    outTruncated.insert(  outTruncated.end(),  outoutTruncated.begin(),  outoutTruncated.end()  );
    outTruncated.insert(  outTruncated.end(),  outoutTruncatedRC.begin(),  outoutTruncatedRC.end()  );
    
  } // end loop over Genome Sequences
  
  
  // Return av data
  Rcpp::DataFrame outDF;
  outDF =
  Rcpp::DataFrame::create(
    Rcpp::Named("Seqid")      = outSeqid,
    Rcpp::Named("Strand")     = outStrand,
    Rcpp::Named("Start")      = outStart,
    Rcpp::Named("End")        = outEnd,
    Rcpp::Named("Truncated")  = outTruncated);
    
  return outDF;
}
