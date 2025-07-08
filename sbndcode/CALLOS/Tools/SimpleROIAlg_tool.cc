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
#include <TROOT.h>
// #include <TVirtualFFT.h>
#include <TMath.h>

namespace callos {

class SimpleROIAlg : public ROIFINDERALG {
public:
  explicit SimpleROIAlg(fhicl::ParameterSet const& p);

  ~SimpleROIAlg() {}

  // Required functions.
  bool ProcessWaveform(std::vector<float> & wvf, std::vector<SimpleROI> & ROI, const int wfChannel);

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
  // std::vector<float> calculateFFT(std::vector<float> &wvf);
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
  int pretrigger = 0.17*wvf.size();
  int segment_size = pretrigger/fBaselineEstimationSegments;
  
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

// std::vector<float> callos::SimpleROIAlg::calculateFFT(std::vector<float> &wvf) {
//   // Calcular la FFT de la waveform
//   int N = wvf.size();

//   // Crear un arreglo para la entrada
//   TVirtualFFT *fft = TVirtualFFT::FFT(1, &N, "R2C EX K");

//   // Copiar los datos del vector de entrada al arreglo de entrada
//   for (int i = 0; i < N; ++i) {
//       fft->SetPoint(i, wvf[i]);
//   }

//   // Ejecutar la FFT
//   fft->Transform();

//   // Obtener los resultados
//   std::vector<float> realPart(N/2 + 1);
//   double re, im;
//   for (int i = 0; i < N/2 + 1; ++i) {
//       fft->GetPointComplex(i, re, im);
//       realPart[i] = re;
//   }

//   // Limpiar
//   delete fft;

//   return realPart;
// }
bool callos::SimpleROIAlg::ProcessWaveform(std::vector<float> & wvf ,std::vector<SimpleROI> & ROI, const int wfChannel)
{
  if (fDefaultROI)
  {

    int wvf_size=wvf.size();
    // Find the best baseline and std from the first 1/10 of the waveform
    auto [best_baseline, best_std] = FindBestBaselineAndSTD(wvf);

    // fConstantBaseline = EstimateBaseline(wave_baseline);
    // if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: baseline is "<<best_baseline<<std::endl;
    // // double baseline_std = BaselineSTD(wave_baseline, fConstantBaseline);
    // if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: baseline std is "<<best_std<<std::endl;
    SubtractBaseline(wvf,best_baseline);
    auto wvf_max = *std::max_element(wvf.begin(), wvf.end());
    auto wvf_min = *std::min_element(wvf.begin(), wvf.end());
    //find the position of the max of the waveform
    // auto max_pos = std::distance(wvf.begin(), std::max_element(wvf.begin(), wvf.end()));
    // if (max_pos > 0.5*wvf_size) return false;
    //auto full_wvf_max = *std::max_element(wvf.begin(), wvf.end());

    // if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: max of the waveform is "<<wvf_max<<std::endl;
    // if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: min of the waveform is "<<wvf_min<<std::endl;
    int readout_window = 768;
    //check max of the waveform is below ThresholdHigh and above ThresholdLow, else there is no signal
    if ((wvf_max > fThresholdHigh ) || (wvf_max<fThresholdLow) || (wvf_min < fThresholdUndershoot) || (wvf_size < readout_window )) return false;

    //If find a noise only window save its charge
    // if (wvf_max < 2 && wvf_min > -2) {
    //   if ((200 + fFinalWindow) > wvf_size) return false;
    //   double aux_charge = std::accumulate(wvf.begin()+200, wvf.begin()+ 200 + fFinalWindow, 0.0);
    //   SimpleROI roi (2*wvf_size, 1);
    //   roi.SetStartTick(0);
    //   roi.SetBaseline(best_baseline);
    //   roi.SetBaselineSTD(best_std);
    //   roi.SetPeak(0);
    //   roi.SetPeakTick(0);
    //   roi.SetEndTick(0);
    //   roi.SetWaveform(wvf);
    //   roi.SetCharge(aux_charge);
    //   // roi.SetFFT_ROI(roi_fft);
    //   ROI.push_back(roi);
    //   return true;
    // }

    // Default ROI finder
    // Signal could be empty or contain single PEs. Read the wvf, look for  peaks above a threshold,
    // compute the charge untill the waveform is below fConstantBaseline again to left and right
    // store the ROI in the ROI vector
    // if(fDebug)
    // {
    //  std::cout << "Baseline: "<<best_baseline<<std::endl;
    // }
    if (fPeakWindow)
    {
      double pretrigger = 0.16*wvf.size(); // 20% of the waveform is the pretrigger
      double pretrigger_max = *std::max_element(wvf.begin(), wvf.begin()+pretrigger);
      double pretrigger_min = *std::min_element(wvf.begin(), wvf.begin()+pretrigger);
      if ((pretrigger_max > 2) || (pretrigger_min < -2)) return false;
      

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
                  if (((wvf[i-k] - std::ceil(best_std)) <= 0))
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
              while ( ((wvf[j] - std::ceil(best_std)) > 0) && (j< (int) wvf.size()) )
              {

                j++;
              }
              int last_bin = j;
              int width = last_bin - first_bin;
              // Width condition (to avoid single spikes) and waveform return baseline in the readout window
              if (width < 15 ){
                i = last_bin;
                
                if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: ROI too small."<< " Width:"<< width <<std::endl;
                continue;
                // return false;
              }
              
              double roi_max = *std::max_element(wvf.begin()+first_bin, wvf.begin()+last_bin);
              if ((roi_max > fRoi_max ) || (roi_max < fRoi_min)){
                i = last_bin;
                continue;
                // return false;
              }
              auto max_pos = std::distance(wvf.begin(), std::max_element(wvf.begin()+first_bin, wvf.begin()+last_bin));
              SimpleROI roi (2*wvf.size(), 1);





              //Align --> Compute the max of the derivate of the rising edge, then do a linear fit near this point and select where it crosses the baseline
              std::vector<float> slope;
              int bins_after_first_max = 10;
              // int derivate_step = 5;
              
              for (int k = 0; k < 30; ++k) {
                  // slope.push_back(wvf[max_pos - k + bins_after_first_max] - wvf[max_pos - k - derivate_step + bins_after_first_max]);
                  slope.push_back((1/12)*(wvf[max_pos - k + bins_after_first_max -2]- 8*wvf[max_pos - k + bins_after_first_max - 1] + 8*wvf[max_pos - k + bins_after_first_max + 1]- wvf[max_pos - k + bins_after_first_max +2]));
              }
              double max_slope = *std::max_element(slope.begin(), slope.end());
              if (max_slope > 4000){
                i = last_bin;
                if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: slope is too big"<<std::endl;
                continue;
                // return false;
              }
              int align = std::distance(slope.begin(), std::max_element(slope.begin(), slope.end()));
              // int align_pos = max_pos - align + bins_after_first_max - std::ceil(derivate_step/2);
              int align_pos = max_pos - align + bins_after_first_max;
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
              if (m == 0) {
                i = last_bin;
                continue;
                // return false;
              }
              float b = (sum_y - m * sum_x) / n;
              float align_fit_pos = - b / m;

              if (fDebug) {
                  // std::cout << "X_0: " << bins_to_fit_x[0] << std::endl;
                  // std::cout << "Y_0: " << bins_to_fit_y[0] << std::endl;
                  // std::cout << "X_n: " << bins_to_fit_x[n-1] << std::endl;
                  // std::cout << "Y_n: " << bins_to_fit_y[n-1] << std::endl;
                  std::cout << "Slope fit: " << m << std::endl;
                  std::cout << "Intercept fit: " << b << std::endl;
                  std::cout << "Align fit pos: " << align_fit_pos << std::endl;
              }

              int start_roi = align_fit_pos - fInitialWindow;
              int end_roi = align_fit_pos + fFinalWindow;
              if (end_roi > (int) wvf.size()){
                i = last_bin;
                if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: end_roi is out of range"<<std::endl;
                return false;
                // continue; 
              }

              bool wvfm_recover_baseline_before = false;
              //int number_of_bins_in_baseline = 0;
              //check if the waveform recovers the baseline in the 30 bins after the peak
              //Dejarlo como antes mirar detras del pico al principio y al final de la ventana
              int bins_in_baseline_before = 0;
              //añadir aprox al entero mayor std_baseline 
              double baseline_std = std::ceil(best_std);
              for (size_t k = 0; k < 10; k++)
              {
                if (((abs(wvf[start_roi + k]) - baseline_std) <= 0))
                {
                  bins_in_baseline_before++;
                }
              }
              if (bins_in_baseline_before > 6) wvfm_recover_baseline_before = true;
              if (wvfm_recover_baseline_before == false) {
                i = last_bin;
                if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: waveform does not recover baseline before the peak"<<std::endl;
                continue;
                // return false;
                
              }
              //double max_full_roi = *std::max_element(wvf.begin()+start_roi, wvf.begin()+end_roi);

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
              if (wvfm_recover_baseline_after == false) {
                i = last_bin;
                if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: waveform does not recover baseline after the peak"<<std::endl;
                continue;
                // return false;
                
              }
              // Lets check the width of the peak
              //Symplify when the wvfm touch the std bejore and after the peak
              // int first_bin_baseline = 0;
              // int last_bin_baseline = 0;
              // for (size_t k = 0; k < 80; k++)
              // {
              //   if (((abs(wvf[max_pos - k]) - baseline_std) <= 0))
              //   {
              //     last_bin_baseline = k;
              //     break;
              //   }
              // }
              // for (size_t k = 0; k < 200; k++)
              // {
              //   if (((abs(wvf[max_pos + k]) - baseline_std) <= 0) && max_pos + k < wvf.size())
              //   {
              //     first_bin_baseline = k;
              //     break;
              //   }
              // }
              // if (last_bin_baseline + first_bin_baseline < 10) return false;

              //Compute the mean peak amplitude wich is located 15 (10) bins after the align_fit_pos for HS (VE) channels by averaging the 2 bins around the peak
              std::vector<int> HS_XA_Channels = {713, 715, 717, 719, 733, 735, 761, 763}; 

              double peak_amplitude_roi = 0;
              int peak_pos;
              int bins_used_peak = 0;

              if (std::find(HS_XA_Channels.begin(), HS_XA_Channels.end(), wfChannel) != HS_XA_Channels.end()) {
                  peak_pos = std::round(align_fit_pos) + 15;
              } else {
                  peak_pos = std::round(align_fit_pos) + 10;
              }

              for (int k = peak_pos - 2; k < peak_pos + 3; k++) {
                  peak_amplitude_roi += wvf[k];
                  bins_used_peak++;
              }

              if (bins_used_peak > 0) {
                  peak_amplitude_roi /= bins_used_peak;
              }


              float aux_charge = 0;

              // roi.SetSmoothROI(smoothed_wvf);
              // roi.SetDerivateROI(roi_derivative);


              for (int k = align_fit_pos - 5; k < align_fit_pos + fFinalWindow; k++) 
              {

                // Add the bins after the threshold until the waveform crosses the baseline
                aux_charge += wvf[k];
                
              }

              if (abs(aux_charge) > 1e6) //Sanity cut to avoid bad charge computation
              {
                i = last_bin;
                if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: aux_charge is too high"<<std::endl;
                continue;
                // return false;
              }
                // if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: found a peak at "<<i<<" with charge "<<aux_charge<<std::endl;
              // int distance_to_end = wvf.size() - end_roi;


              double peak_charge = 0;
              for (int k = align_fit_pos; k < align_fit_pos + 50; k++) //De normal +30
              {
                peak_charge += wvf[k];
              }
              if (peak_charge < 0.70*aux_charge){
                i = last_bin;
                if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: peak charge is too spread"<<std::endl;
                continue;
                // return false;
                
              }

            

              std::vector<float> roi_wvf; // Declare roi_wvf here
              // if (distance_to_end < 30)
              // { 
              roi_wvf.assign(wvf.begin() + start_roi - 10, wvf.begin() + end_roi + 10);
                //roi_wvf.assign(wvf.begin() + align_fit_pos - 20, wvf.begin() + align_fit_pos + 130);
              // }
             
              // else
              // {
              //   roi_wvf.assign(wvf.begin() + start_roi - 8, wvf.begin() + end_roi + 30);
              //   //roi_wvf.assign(wvf.begin() + align_fit_pos - 20, wvf.begin() + align_fit_pos + 130);
              // }
              //Condition to get the mean ampltitudo of the 1 PE peak
              //if (aux_charge < 35 || aux_charge > 70) return false;
              //Check if ROI waveform is below roi_max
              for (size_t k = 0; k < roi_wvf.size(); k++)
              {
                if (abs(roi_wvf[k]) > fRoi_max){
                i = last_bin;
                continue;
                // return false;
                }
              }

              //Compute FFT of the ROI
              // std::vector<float> roi_fft = calculateFFT(wvf);

              // SubtractBaseline(roi_wvf, best_baseline);

              roi.SetStartTick(start_roi);
              roi.SetBaseline(best_baseline);
              roi.SetBaselineSTD(best_std);
              roi.SetPeak(peak_amplitude_roi);
              roi.SetPeakTick(align_fit_pos);
              roi.SetEndTick(end_roi);
              roi.SetWaveform(roi_wvf);
              roi.SetCharge(aux_charge);
              // roi.SetFFT_ROI(roi_fft);
              ROI.push_back(roi);
              i = end_roi;
            }
          else
          {
            if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform: waveform does not recover baseline before the peak"<<std::endl;
            continue;
            // return false;
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
    SimpleROI roi (wvf.size(), 1);
    // double peak_charge = 0;
    // double window_charge = 0;
    // for (int k = 0; k < 768; k++) 
    // {
    //   window_charge += wvf[k];
    // }
    // for (int k = 100; k < 300; k++) 
    // {
    //   peak_charge += wvf[k];
    // }
    // if (peak_charge < 0.50*window_charge) return false;
    double pretrigger = 0.1*wvf.size(); // 10% of the waveform is the pretrigger
    for (size_t i = pretrigger; i < wvf.size(); i++) 
    {
      
      bool wvfm_recover_baseline = false;
      int bins_before = 0;
      //int number_of_bins_in_baseline = 0;
      if (wvf[i] > fThresholdLow ) //found a peak?
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
            
            //double roi_max = *std::max_element(wvf.begin()+first_bin, wvf.begin()+last_bin);
            auto max_pos = std::distance(wvf.begin(), std::max_element(wvf.begin()+first_bin, wvf.begin()+last_bin));





            //Align --> Compute the max of the derivate of the rising edge, then do a linear fit near tis point and select where it crosses the baseline
            std::vector<float> slope;
            int bins_after_first_max = 10;
            int derivate_step = 5;
            for (int k = 0; k < 50; ++k) {
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
            for (int k = align_pos - 5; k < align_pos + 5 ; ++k) {
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
            if (m == 0) return false;
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

            int start_roi = align_fit_pos - fInitialWindow;
            int end_roi = align_fit_pos + fFinalWindow;
            if (end_roi > (int) wvf.size()) return false;
            bool wvfm_recover_baseline_before = false;

            int bins_in_baseline_before = 0;
            //añadir aprox a l entero mayor std_baseline 
            for (size_t k = 0; k < 100; k++)
            {
              if (((abs(wvf[start_roi + k])) <= 1))
              {
                bins_in_baseline_before++;
              }
            }
            if (bins_in_baseline_before > 66) wvfm_recover_baseline_before = true;
            if (wvfm_recover_baseline_before == false) return false;
            bool wvfm_recover_baseline_after = false;
            // //check if the waveform recovers the baseline in the 30 bins after the peak
            int bins_in_baseline_after = 0;
            // if (end_roi + 10 > (int) wvf.size()) return false;
            for (size_t k = 0; k < 10; k++)
            {
              if (((abs(wvf[end_roi - k])) <= 1))
              {
                bins_in_baseline_after++;
              }
            }
            if (bins_in_baseline_after > 6) wvfm_recover_baseline_after = true;
            if (wvfm_recover_baseline_after == false) return false;
            std::vector<float> roi_wvf; // Declare roi_wvf here
              // if (distance_to_end < 30)
              // { 
            roi_wvf.assign(wvf.begin() + start_roi, wvf.begin() + end_roi);
            i = end_roi;
            roi.SetWaveform(roi_wvf);
            roi.SetBaseline(best_baseline);
            roi.SetPeak(wvf_max);
            ROI.push_back(roi);
          }
        }
    }




 
  }   
  
   //if (fDebug) std::cout << "SimpleROIAlg::ProcessWaveform"<<std::endl;
  return true;

}

DEFINE_ART_CLASS_TOOL(callos::SimpleROIAlg)

