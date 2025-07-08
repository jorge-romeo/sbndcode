////////////////////////////////////////////////////////////////////////
// Specific class tool for ROIFindetAlg tool
// File: SimpleROI_tool.cc
// Base class:        ROIFINDERALG.hh
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

#include "sbndcode/CALLOS/Tools/ROIFinderAlg.hh"
#include <iostream>
#include <vector>
#include <numeric>  // Incluir este encabezado para std::accumulate
#include <cmath>
#include <algorithm>

namespace callos {

class SimpleROIAlg : public ROIFINDERALG {
public:
  explicit SimpleROIAlg(fhicl::ParameterSet const& p);

  ~SimpleROIAlg() {}

  // Required functions.
  bool ProcessWaveform(std::vector<float> & wvf, std::vector<SimpleROI> & ROI);

private:
  // Declare member data here.
  bool fDebug;
  float fConstantBaseline;
  float fThresholdLow;
  float fThresholdHigh;
  float fThresholdUndershoot;
  int fBinsBefore;
  int fMaxADC;
  bool fDefaultROI;
  int fInitialWindow;
  int fFinalWindow;
  int fBaselineEstimationSegments;
  bool fPeakWindow;
  double fRoi_max;
  double fRoi_min;

  // Declare member functions
  void SubtractBaseline(std::vector<float> &wf, float baseline);
  void SubtractBaseline(std::vector<float> &wf, std::vector<float> &baseline);  // overload for non-constant baselines
  float EstimateBaseline(std::vector<float> &wf);  // Cambiado a float
  float BaselineSTD(std::vector<float> &wf, float baseline);  // Cambiado a float
  std::pair<float, float> FindBestBaselineAndSTD(const std::vector<float> &wvf);
};

} // namespace callos


callos::SimpleROIAlg::SimpleROIAlg(fhicl::ParameterSet const& p)
{
  fDebug = p.get< bool >("Debug",false);
  fConstantBaseline = p.get< float >( "ConstantBaseline", 0.0);
  fThresholdLow = p.get< float >( "ThresholdLow", 0.0);
  fThresholdHigh = p.get< float >( "ThresholdHigh", 0.0);
  fBinsBefore = p.get< int >( "BinsBefore", 0);
  fThresholdUndershoot = p.get< float >( "ThresholdUndershoot", -4);
  fMaxADC = p.get< float >( "MaxADC", 4000);
  fDefaultROI = p.get< bool >( "DefaultROI", true);
  fInitialWindow = p.get< int >( "InitialWindow", 0);
  fFinalWindow = p.get< int >( "FinalWindow", 0);
  fBaselineEstimationSegments = p.get< int >( "BaselineEstimationSegments", 2);
  fPeakWindow = p.get< bool >( "PeakWindow", false);
  fRoi_max = p.get< double >( "Roi_max", 0);
  fRoi_min = p.get< double >( "Roi_min", 0);

  if (fDebug) 
  {
    std::cout << "SimpleROIAlg fhicl congif::"<<std::endl;
    std::cout << "ConstantBaseline: "<<fConstantBaseline<<std::endl;
    std::cout << "ThresholdLow: "<<fThresholdLow<<std::endl;
    std::cout << "ThresholdHigh: "<<fThresholdHigh<<std::endl;
  }

}

// constant baseline substractor
void callos::SimpleROIAlg::SubtractBaseline(std::vector<float> &wf, float baseline){
  std::transform (wf.begin(), wf.end(), wf.begin(), [&](float _x) { return _x-baseline; }  );
  return;
}

// compute the baseline of the waveform
float callos::SimpleROIAlg::EstimateBaseline(std::vector<float> &wf) {
  return std::accumulate(wf.begin(), wf.end(), 0.0f) / wf.size();
}

// compute the standard deviation of the baseline
float callos::SimpleROIAlg::BaselineSTD(std::vector<float> &wf, float baseline) {
  double mean = baseline;
  double variance = 0.0;
  for (const auto& value : wf) {
    variance += (value - mean) * (value - mean);
  }

  if (wf.size() > 0)
  { 
    variance /= wf.size();
  }
  else
  {
    // std::cout << "SimpleROIAlg::BaselineSTD: wf.size() is 0"<<std::endl;
    return 0.0;
  }
  return std::sqrt(variance);
}

// Function to find the baseline and std of the segment with the lowest std
std::pair<float, float> callos::SimpleROIAlg::FindBestBaselineAndSTD(const std::vector<float> &wvf) {
  int tenth = wvf.size() / 10;
  int segment_size = tenth/fBaselineEstimationSegments;
  
  float best_baseline = 0.0f;
  float best_std = std::numeric_limits<float>::max();

  for (int i = 0; i < fBaselineEstimationSegments; ++i) {
    std::vector<float> segment(wvf.begin() + i * segment_size, wvf.begin() + (i + 1) * segment_size);
    float baseline = EstimateBaseline(segment);
    float std = BaselineSTD(segment, baseline);
    // std::cout << "SimpleROIAlg::ProcessWaveform: intermediate baseline is "<<baseline<<std::endl;
    // std::cout << "SimpleROIAlg::ProcessWaveform: intermediate baseline std is "<<std<<std::endl;
    if (std < best_std) {
      best_std = std;
      best_baseline = baseline;
    }
  }

  return {best_baseline, best_std};
}
// // variable baseline subtractor 
// void callos::SimpleROIAlg::SubtractBaseline(std::vector<float> &wf, std::vector<float> &baseline){
//   std::transform (wf.begin(), wf.end(), baseline.begin(), wf.begin(), std::minus<float>());
//   return;
// }

bool callos::SimpleROIAlg::ProcessWaveform(std::vector<float> & wvf ,std::vector<SimpleROI> & ROI)
{
  if (fDefaultROI)
  {
    // constant for now, will update 


    int wvf_size=wvf.size();
    // Find the best baseline and std from the first 1/10 of the waveform
    auto [best_baseline, best_std] = FindBestBaselineAndSTD(wvf);

    // fConstantBaseline = EstimateBaseline(wave_baseline);
    if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: baseline is "<<best_baseline<<std::endl;
    // double baseline_std = BaselineSTD(wave_baseline, fConstantBaseline);
    if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: baseline std is "<<best_std<<std::endl;
    SubtractBaseline(wvf,best_baseline);
    auto wvf_max = *std::max_element(wvf.begin(), wvf.end());
    auto wvf_min = *std::min_element(wvf.begin(), wvf.end());
    //find the position of the max of the waveform
    // auto max_pos = std::distance(wvf.begin(), std::max_element(wvf.begin(), wvf.end()));
    // if (max_pos > 0.5*wvf_size) return false;
    //auto full_wvf_max = *std::max_element(wvf.begin(), wvf.end());

    if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: max of the waveform is "<<wvf_max<<std::endl;
    if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: min of the waveform is "<<wvf_min<<std::endl;
    int readout_window = 768;
    //check max of the waveform is below ThresholdHigh and above ThresholdLow, else there is no signal
    if ((wvf_max > fThresholdHigh ) || (wvf_max<fThresholdLow) || (wvf_min < fThresholdUndershoot) || (wvf_size < readout_window )) return false;

    // Default ROI finder
    // Signal could be empty or contain single PEs. Read the wvf, look for  peaks above a threshold,
    // compute the charge untill the waveform is below fConstantBaseline again to left and right
    // store the ROI in the ROI vector
    if(fDebug)
    {
     std::cout << "Baseline: "<<best_baseline<<std::endl;
    }
    if (fPeakWindow)
    {
      double pretrigger = 0.1*wvf.size(); // 10% of the waveform is the pretrigger
      for (size_t i = pretrigger; i < wvf.size(); i++) 
      {
        
        bool wvfm_recover_baseline = false;
        int bins_before = 0;
        //int number_of_bins_in_baseline = 0;
        if (wvf[i] > fRoi_min ) //found a peak?
        { 
          // condition to find a peak, 50 bins before the waveform should be equal or below the constant baseline
          //ind the first bin before the peak that is near the constant baseline and set it as the first bin of the ROI
            for (size_t k = 0; k < 50; k++)
            {
                  if (((wvf[i-k] - best_std) <= 0))
                  {
                    wvfm_recover_baseline = true;
                    bins_before = k;
                    break;
                  }
            }
          
          // Compute the first bin of the ROI
            if (wvfm_recover_baseline) 
            {
              int first_bin = i-bins_before;
              int j = i;
              // Add the bins after the threshold until the waveform crosses the baseline
              while ( ((wvf[j] - best_std) > 0) && (j< (int) wvf.size()) )
              {

                j++;
              }
              int last_bin = j;
              // Width condition (to avoid single spikes) and waveform return baseline in the readout window
              if (last_bin - first_bin < 10 ) return false;
              i = j;
              double roi_max = *std::max_element(wvf.begin()+first_bin, wvf.begin()+last_bin);
              if ((roi_max > fRoi_max ) || (roi_max < fRoi_min)) return false;



              auto max_pos = std::distance(wvf.begin(), std::max_element(wvf.begin()+first_bin, wvf.begin()+last_bin));
              SimpleROI roi (wvf.size(), 1);
              int start_roi = max_pos - fInitialWindow;
              int end_roi = max_pos + fFinalWindow;
              if (end_roi > (int) wvf.size()) return false;



              bool wvfm_recover_baseline_before = false;
              //int number_of_bins_in_baseline = 0;
              //check if the waveform recovers the baseline in the 30 bins after the peak
              //Dejarlo como antes mirar detras del pico al principio y al final de la ventana
              int bins_in_baseline_before = 0;
              //añadir aprox a l entero mayor std_baseline 
              double baseline_std = std::ceil(best_std);
              for (size_t k = 0; k < 10; k++)
              {
                if (((abs(wvf[start_roi + k]) - baseline_std) <= 0))
                {
                  bins_in_baseline_before++;
                }
              }
              if (bins_in_baseline_before > 6) wvfm_recover_baseline_before = true;
              if (wvfm_recover_baseline_before == false) return false;
              double max_full_roi = *std::max_element(wvf.begin()+start_roi, wvf.begin()+end_roi);

              bool wvfm_recover_baseline_after = false;
              // //check if the waveform recovers the baseline in the 30 bins after the peak
              int bins_in_baseline_after = 0;
              // if (end_roi + 10 > (int) wvf.size()) return false;
              for (size_t k = 0; k < 10; k++)
              {
                if (((abs(wvf[end_roi - k]) - baseline_std) <= 0))
                {
                  bins_in_baseline_after++;
                }
              }
              if (bins_in_baseline_after > 6) wvfm_recover_baseline_after = true;
              if (wvfm_recover_baseline_after == false) return false;

              // Lets check the width of the peak
              //Symplify when the wvfm touch the std bejore and after the peak
              int first_bin_baseline = 0;
              int last_bin_baseline = 0;
              for (size_t k = 0; k < 80; k++)
              {
                if (((abs(wvf[max_pos - k]) - baseline_std) <= 0))
                {
                  last_bin_baseline = k;
                  break;
                }
              }
              for (size_t k = 0; k < 200; k++)
              {
                if (((abs(wvf[max_pos + k]) - baseline_std) <= 0) && max_pos + k < wvf.size())
                {
                  first_bin_baseline = k;
                  break;
                }
              }
              if (last_bin_baseline + first_bin_baseline < 10) return false;

              roi.SetCharge(0);
              roi.SetStartTick(start_roi);
              roi.SetBaseline(best_baseline);
              roi.SetBaselineSTD(best_std);
              roi.SetPeak(max_full_roi);
              roi.SetPeakTick(max_pos);
              float aux_charge = 0;
              roi.SetEndTick(end_roi);
              // roi.SetSmoothROI(smoothed_wvf);
              // roi.SetDerivateROI(roi_derivative);


              for (int k = start_roi; k < end_roi; k++) 
              {

                // Add the bins after the threshold until the waveform crosses the baseline
                aux_charge += wvf[k];
                
              }
                // if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: found a peak at "<<i<<" with charge "<<aux_charge<<std::endl;
              // int distance_to_end = wvf.size() - end_roi;


              double peak_charge = 0;
              for (int k = max_pos - 15; k < max_pos + 15; k++) 
              {
                peak_charge += wvf[k];
              }
              if (peak_charge < 0.50*aux_charge) return false;
              
              //Align --> Compute the max of the derivate of the rising edge, then do a linear fit near tis point and select where it crosses the baseline
              std::vector<float> slope;
              int bins_after_first_max = 20;
              int derivate_step = 5;
              for (int k = 0; k < 40; ++k) {
                  slope.push_back(wvf[max_pos - k + bins_after_first_max] - wvf[max_pos - k - derivate_step + bins_after_first_max]);
              }
              int align = std::distance(slope.begin(), std::max_element(slope.begin(), slope.end()));
              int align_pos = max_pos - align + bins_after_first_max - std::ceil(derivate_step/2);
              if (fDebug) {
                  std::cout << "Slope: ";
                  for (const auto& val : slope) {
                      std::cout << val << " ";
                  }
                  std::cout << std::endl;
                  std::cout << "Align: " << align << std::endl;
                  std::cout << "Align pos: " << align_pos  << std::endl;
              }
              
              //Linear fit to linear edge
              std::vector<float> bins_to_fit_x;
              std::vector<float> bins_to_fit_y;
              for (int k = align_pos - 2; k < align_pos + 3 ; ++k) {
                  bins_to_fit_x.push_back(k);
                  bins_to_fit_y.push_back(wvf[k]);
              }
              size_t n = bins_to_fit_x.size();
              double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0;

              for (size_t k = 0; k < n; ++k) {
                  sum_x += bins_to_fit_x[k];
                  sum_y += bins_to_fit_y[k];
                  sum_xy += bins_to_fit_x[k] * bins_to_fit_y[k];
                  sum_x2 += bins_to_fit_x[k] * bins_to_fit_x[k];
              }

              // Calcular la pendiente (m) y el intercepto (b)
              float m = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x * sum_x);
              float b = (sum_y - m * sum_x) / n;
              float align_fit_pos = - b / m;
              if (fDebug) {
                  std::cout << "X_0: " << bins_to_fit_x[0] << std::endl;
                  std::cout << "Y_0: " << bins_to_fit_y[0] << std::endl;
                  std::cout << "X_n: " << bins_to_fit_x[n-1] << std::endl;
                  std::cout << "Y_n: " << bins_to_fit_y[n-1] << std::endl;
                  std::cout << "Slope fit: " << m << std::endl;
                  std::cout << "Intercept fit: " << b << std::endl;
                  std::cout << "Align fit pos: " << align_fit_pos << std::endl;
              }


            

              std::vector<float> roi_wvf; // Declare roi_wvf here
              // if (distance_to_end < 30)
              // { 
                //roi_wvf.assign(wvf.begin() + start_roi - 20, wvf.begin() + end_roi);
              roi_wvf.assign(wvf.begin() + align_fit_pos - 20, wvf.begin() + align_fit_pos + 160);
              // }
             
              // else
              // {
              //   //roi_wvf.assign(wvf.begin() + start_roi - 20, wvf.begin() + end_roi + 30);
              //   roi_wvf.assign(wvf.begin() + align_fit_pos - 20, wvf.begin() + align_fit_pos + 130);
              // }
              //Condition to get the mean ampltitudo of the 1 PE peak
              //if (aux_charge < 35 || aux_charge > 70) return false;
              //Check if ROI waveform is below roi_max
              for (size_t k = 0; k < roi_wvf.size(); k++)
              {
                if (abs(roi_wvf[k]) > fRoi_max) return false;
              }
              // SubtractBaseline(roi_wvf, best_baseline);
              roi.SetWaveform(roi_wvf);
              roi.SetCharge(aux_charge);
              ROI.push_back(roi);
              i = end_roi;
            }     
        }
      }  

      
    }
    else{
      SimpleROI roi (wvf.size(), 1);
      //this should be updated with something more meaningful later on
      roi.SetCharge(0);
      roi.SetStartTick(fInitialWindow);
      roi.SetBaseline(best_baseline);
      roi.SetBaselineSTD(best_std);
      float aux_charge = 0;
      roi.SetEndTick(fFinalWindow);
      //roi.SetAverageWaveform(wvf);
      for (int i = fInitialWindow; i < fFinalWindow; i++) 
      {

        // Add the bins after the threshold until the waveform crosses the baseline
        aux_charge += wvf[i];
        
      }
        // if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: found a peak at "<<i<<" with charge "<<aux_charge<<std::endl;
      std::vector<float> roi_wvf(wvf.begin()+ fInitialWindow, wvf.begin() + fFinalWindow);
      // SubtractBaseline(roi_wvf, best_baseline);
      roi.SetWaveform(roi_wvf);
      roi.SetCharge(aux_charge);
      ROI.push_back(roi);
    }

  }
  else{
    // constant for now, will update 
    std::vector<float> wave;
    //std::vector<float> wave_baseline;
    // look at the second half of the waveform (scintillation tail, after the prompt light is gone)
    int wvf_size=wvf.size();
    int half=2*wvf.size()/3;
    wave.reserve(half);
    wave.assign(wvf.begin()+half, wvf.end());


    // Find the best baseline and std from the first 1/10 of the waveform
    auto [best_baseline, best_std] = FindBestBaselineAndSTD(wvf);
    SubtractBaseline(wvf,best_baseline);
    // fConstantBaseline = EstimateBaseline(wave_baseline);
    if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: baseline is "<<best_baseline<<std::endl;
    // double baseline_std = BaselineSTD(wave_baseline, fConstantBaseline);
    if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: baseline std is "<<best_std<<std::endl;
    SubtractBaseline(wave,best_baseline);
    auto wvf_max = *std::max_element(wave.begin(), wave.end());
    auto wvf_min = *std::min_element(wave.begin(), wave.end());
    auto full_wvf_max = *std::max_element(wvf.begin(), wvf.end()) + best_baseline;

    if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: max of the waveform is "<<wvf_max<<std::endl;
    if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: min of the waveform is "<<wvf_min<<std::endl;
    int readout_window = 768;
    //check max of the waveform is below ThresholdHigh and above ThresholdLow, else there is no signal
    if ((wvf_max > fThresholdHigh ) || (wvf_max<fThresholdLow) || (wvf_min < fThresholdUndershoot) || (full_wvf_max >= fMaxADC )|| (wvf_size < readout_window )) return false;

    if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: passed max threshold requisite"<<std::endl;
    // Signal could be empty or contain single PEs. Read the wvf, look for  peaks above a threshold,
    // compute the charge untill the waveform is below fConstantBaseline again to left and right
    // store the ROI in the ROI vector

    for (size_t i = 0; i < wave.size(); i++) {
      
      
      bool wvfm_recover_baseline = false;
      int bins_before = 0;
      //int number_of_bins_in_baseline = 0;
      if (wave[i] > fThresholdLow ) { //found a peak
        // condition to find a peak, 50 bins before the waveform should be equal or below the constant baseline
        //ind the first bin before the peak that is near the constant baseline and set it as the first bin of the ROI
          for (size_t j = 0; j < 50; j++)
          {
                if (((wvf[half + i - j] - best_std) <= 0))
                {
                  // for (size_t k = 0; k < 10; k++) 
                  // {
                  //   if ((wvf[half + i - j - k] - best_std) <= best_baseline) 
                  //   {
                  //     number_of_bins_in_baseline++;
                  //   }
                  //   if (number_of_bins_in_baseline > 4) 
                  //   {
                      wvfm_recover_baseline = true;
                      bins_before = j;
                      break;
                  //   }
                  // }
                }
          }
        
        
        // Compute the first bin of the ROI
          if (wvfm_recover_baseline) 
          {
              int first_bin = i+half-bins_before;
              SimpleROI roi(wave.size(), 1);//this should be updated with something more meaningful later on
              roi.SetCharge(0);
              roi.SetStartTick(first_bin);
              
              roi.SetBaseline(best_baseline);
              roi.SetBaselineSTD(best_std);

              float aux_charge = 0;
              int j = i;
              // Add the bins after the threshold until the waveform crosses the baseline
              while ( ((wave[j] - best_std) > 0) && (j< (int) wave.size()) )
              {
                // roi[j-i] = wave[j];
                aux_charge += wave[j];
                j++;
              }
              int bins_in_baseline_after_roi = 0;
              // for (int z = 0; z < 10; z++) 
              // {
              //   if ((wave[j+z] - best_std) > 0)
              //   {
              //     bins_in_baseline_after_roi++;
              //     aux_charge += wave[j+z];
              //   }
              // }
              int last_bin = j+half+bins_in_baseline_after_roi;
              roi.SetEndTick(last_bin);
              // Add the bins before the threshold
              for (int k = 0; k < bins_before; k++) 
              {
                if ((i-k) < 0) break;

                if (wvf[half+i-k] > best_baseline) 
                {
                  // roi[j-i+k] = wave[i-k];
                  aux_charge += wvf[half+i-k];
                }
              }
              if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: found a peak at "<<i<<" with charge "<<aux_charge<<std::endl;
              // Width condition (to avoid single spikes) and waveform return baseline in the readout window
              if ((last_bin - first_bin > 10) && (last_bin - first_bin < 150) && (wave[j] <= 0) && (j < (int) wave.size())) 
              {
                std::vector<float> roi_wvf(wvf.begin()+ first_bin -20, wvf.begin() + last_bin+ 20);
                //Find the maximum of the ROI requiring at least 3 bins with this value
                double roi_max = *std::max_element(roi_wvf.begin()+20, roi_wvf.end()-20);
                //int roi_max_count = std::count(roi_wvf.begin()+20, roi_wvf.end()-20, roi_max);
                // if (roi_max_count > 3)
                // {
                //    continue;
                // }   
                // else
                // {
                //   roi_max = roi_max - 1;
                // }  
                roi.SetPeak(roi_max);
                roi.SetWaveform(roi_wvf);
                roi.SetCharge(aux_charge);
                ROI.push_back(roi);
              }
              // roi.SetCharge(aux_charge);
              // ROI.push_back(roi);
              i = j;//skip the waveform, find the next peak
          
            }
      }
    }
  }   
  
   //if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform"<<std::endl;
  return true;

}

DEFINE_ART_CLASS_TOOL(callos::SimpleROIAlg)

