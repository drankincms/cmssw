#include "HeterogeneousCore/SonicTriton/interface/TritonConverterBase.h"
#include "HeterogeneousCore/SonicTriton/interface/TritonData.h"

class FloatStandardConverter : public TritonConverterBase<float> {
public:
  FloatStandardConverter(const edm::ParameterSet& conf)
      : TritonConverterBase<float>(conf) {}

  uint8_t* convert(const float* in) const override;
};

DEFINE_EDM_PLUGIN(TritonConverterFactory<float>, FloatStandardConverter, "FloatStandardConverter");

uint8_t* FloatStandardConverter::convert(const float *in) {
    return reinterpret_cast<const uint8_t*>(in);
}
