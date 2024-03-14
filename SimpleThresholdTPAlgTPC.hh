#ifndef DUNEANA_TRIGGERSIM_SIMPLETHRESHOLDTPALGTPC_hh
#define DUNEANA_TRIGGERSIM_SIMPLETHRESHOLDTPALGTPC_hh

#include "fhiclcpp/ParameterSet.h" 

#include "duneana/TriggerSim/TPAlgTools/TPAlgTPCTool.hh"


// Standard (or 'default') TPGen algorithm: pedsub + static threshold hit finding (no filter)
namespace duneana {

  class SimpleThresholdTPAlgTPC : public TPAlgTPCTool {

  public:
    explicit SimpleThresholdTPAlgTPC(fhicl::ParameterSet const& ps) :
      verbosity_(ps.get<int>("verbosity",0)),
      m_contingency(ps.get<short>("contingency",10)),   //contingency limit for pedsub
      m_threshold(ps.get<short>("threshold",40))        //hit finding threshold
    {}


    //pedestal subtraction
    std::vector<short> subtract_pedestal(const std::vector<short> &adcs, 
					 const short m_contingency);
    //hit finding
    void find_hits( const std::vector<short> &adcs, 
		    std::vector<dunedaq::trgdataformats::TriggerPrimitive> &hits);

    //waveform processing --> apply all algorithm parts and output hits 
    void process_waveform(std::vector<short> const& adcs, 
			  dunedaq::trgdataformats::channel_t const channel,
			  dunedaq::trgdataformats::detid_t const detid,
			  dunedaq::trgdataformats::timestamp_t const start_time,
			  std::vector<dunedaq::trgdataformats::TriggerPrimitive> & tps_out); 
        
  private:
    int verbosity_;
    short m_contingency;
    short m_threshold;
  };

  std::vector<short> SimpleThresholdTPAlgTPC::subtract_pedestal( const std::vector<short> &adcs , const short m_contingency)
  {
    short median = adcs[0];                                                           
    std::vector<short> pedestal(adcs.size(), 0);
    short accumulator=0;
    
    for(size_t i=0; i<adcs.size(); ++i){
      short sample = adcs[i]; // current sample                                                                       
    
      if(sample>median) ++accumulator;
      if(sample<median) --accumulator;

      //Update pedestal if contingency limit exceeded by accumulator.
      if(accumulator > m_contingency){
	++median;
	accumulator=0;
      }
      if(accumulator < -1*m_contingency){
	--median;
	accumulator=0;
      }
    
      pedestal[i]= sample - median;
    }
  
    return pedestal;
  }

  void SimpleThresholdTPAlgTPC::find_hits(const std::vector<short> &adcs, std::vector<dunedaq::trgdataformats::TriggerPrimitive> &hits)
  {
    bool is_hit = false;
    bool was_hit = false;
    std::vector<int> hit_charge;

    //initialize the hit 
    dunedaq::trgdataformats::TriggerPrimitive hit;

    
    for(size_t isample=0; isample<adcs.size()-1; ++isample){
      short adc  = adcs[isample];
      is_hit = adc >  (short)m_threshold;
      if(is_hit && !was_hit) {
	hit_charge.push_back(adc); 
	hit.time_start   = isample;
	hit.adc_integral = adc;
	hit.time_over_threshold = 1;
      }
      if(is_hit && was_hit) {
	hit.adc_integral += adc;
	hit.time_over_threshold += 1;
	hit_charge.push_back(adc);
      }
      if(!is_hit && was_hit) {
	//find the peak time and peak ADC
	auto peak_adc_iter = std::max_element(hit_charge.begin(), hit_charge.end());
	hit.adc_peak = *peak_adc_iter; 
	hit.time_peak = std::distance(hit_charge.begin(), peak_adc_iter);

	//Save hit and reset 
	hits.push_back(hit);
	hit_charge.clear();
      }
      was_hit = is_hit;
    }
  }

  void SimpleThresholdTPAlgTPC::process_waveform(std::vector<short> const& adcs, 
					 dunedaq::trgdataformats::channel_t const channel,
					 dunedaq::trgdataformats::detid_t const detid,
					 dunedaq::trgdataformats::timestamp_t const start_time,
					 std::vector<dunedaq::trgdataformats::TriggerPrimitive> & tps_out) 
  {
    // Get pedestal subtracted ADCs
    std::vector<short> pedsub_adcs = subtract_pedestal(adcs, m_contingency);

    // Get hits for this waveform 
    auto hits = std::vector<dunedaq::trgdataformats::TriggerPrimitive>();
    find_hits(pedsub_adcs, hits);
  
    //If signal found - loop over TPs and set global parameters
    if (!hits.empty()) {
      for (size_t j = 0; j< hits.size(); j++){
	dunedaq::trgdataformats::TriggerPrimitive this_tp = hits[j];
	  
	this_tp.channel = channel;
	this_tp.detid = detid;
	
	this_tp.type =      dunedaq::trgdataformats::TriggerPrimitive::Type::kTPC;
	this_tp.algorithm = dunedaq::trgdataformats::TriggerPrimitive::Algorithm::kTPCDefault;
	
	this_tp.flag = 0;
  
	this_tp.time_start = this_tp.time_start + start_time; // is this correct?
	this_tp.time_over_threshold = this_tp.time_over_threshold + start_time; //start_time + adcs.size()*this->ADC_SAMPLING_RATE_IN_DTS;
	this_tp.time_peak = this_tp.time_peak + start_time; //adcs.size()*this->ADC_SAMPLING_RATE_IN_DTS / 2;
      
      
	tps_out.push_back(this_tp);
      }
    }
  }
}

#endif
