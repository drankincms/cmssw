#ifndef __TritonConverterBase_H__
#define __TritonConverterBase_H__

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"

#include <string>

template <typename DT>
class TritonConverterBase {
public:
  TritonConverterBase(const edm::ParameterSet& conf) : _converterName(conf.getParameter<std::string>("converterName")) {}
  TritonConverterBase(const TritonConverterBase&) = delete;
  virtual ~TritonConverterBase() = default;
  TritonConverterBase& operator=(const TritonConverterBase&) = delete;

  virtual uint8_t* convert(const DT* in) const = 0;

  const std::string& name() const { return _converterName; }

private:
  const std::string _converterName;
};

#include "FWCore/PluginManager/interface/PluginFactory.h"
template <typename DT>
using TritonConverterFactory = edmplugin::PluginFactory<TritonConverterBase<DT>*(const edm::ParameterSet&)>;

#endif
