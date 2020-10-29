#include "HeterogeneousCore/SonicTriton/interface/TritonConverterBase.h"

class FloatStandardConverter : public TritonConverterBase<float> {
public:
  FloatStandardConverter(const edm::ParameterSet& conf)
      : TritonConverterBase<float>(conf) {}

  const uint8_t* convert(const float* in) override;
};

DEFINE_EDM_PLUGIN(TritonConverterFactory<float>, FloatStandardConverter, "FloatStandardConverter");

const uint8_t* FloatStandardConverter::convert(const float *in) {
    return reinterpret_cast<const uint8_t*>(in);
}
