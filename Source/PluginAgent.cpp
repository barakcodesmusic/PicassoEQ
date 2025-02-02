#include "PluginAgent.h"

#include "ollama.hpp"

namespace pagent {

PluginAgent::PluginAgent(PluginType&& pluginType) :
	m_pluginType(std::move(pluginType)) 
{
};

PluginAgent::~PluginAgent() {};

void 
PluginAgent::registerParameter(ParameterInfo&& p) 
{
	m_parameterInfos.push_back(std::move(p));
}

std::string 
PluginAgent::inferParameters(std::string& userPrompt) 
{
	ollama::response infer = ollama::generate("deepseek-r1:7b", "Why is the sky blue?");
	return infer.as_simple_string();
}

}

