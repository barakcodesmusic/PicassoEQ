#pragma once

#include <string>
#include <vector>
#include <map>

namespace pagent {

struct ParameterRange {
	float min;
	float max;
};

enum class PluginType {
	EQ,
	Compressor,
	Reverb,
	Delay,
	Chorus,
	Flanger,
	Phaser,
	Filter,
	Other
};

using ParameterName = std::string;
using ParameterValue = float;

struct ParameterInfo {
	ParameterName pname;
	ParameterRange prange;
};

using ParamaterMap = std::map<ParameterName, ParameterValue>;

static std::string MODEL_PROMPT =
	"You are a professional music producer working on your next biggest hit. "
	"You are currently working on a new song, specifically you are tweaking some parameters for a $1 plugin at the moment. "
	"The $1 plugin offers the following parameters in the format <parameter>:<parameter range> for manipulating the input sound: $2"
	"Specifically, you are going for a sound described as: $3."
	"Select parameter values within the given parameter ranges for each given parameter in order to best achieve the sound you have described. "
	"Output your decision in the following format: <parameter1 name>:<parameter1 value>, <parameter2 name>:<parameter2 value>, etc."
	"For instance, if you were working on a reverb with the following parameters: Size:5-50, Decay:0.5ms-60s you might output Size:10, Decay:20s"
	"Importantly, don't include any extra output, just the <parameter name>:<parameter value> list as given in the example. "
	"This output will be fed into another system which expects the input to be formatted as specified in the previous line.";

class PluginAgent {
public:
	PluginAgent(PluginType&& pluginType);
	~PluginAgent();
	void registerParameter(ParameterInfo&& p);
	std::string inferParameters(std::string& userPrompt);

private:
	PluginType m_pluginType;
	std::vector<ParameterInfo> m_parameterInfos;
};

}

