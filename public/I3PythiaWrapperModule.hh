#ifndef I3PYTHIAWRAPPERMODULE_HH
#define I3PYTHIAWRAPPERMODULE_HH

#include "icetray/I3ConditionalModule.h"
#include "icetray/I3Module.h"
#include "Pythia8/Pythia.h"

/// @brief A module that encapsulates Pythia that can be used as an event generator.
class I3PythiaWrapperModule : public I3ConditionalModule
{
public:
    I3PythiaWrapperModule(const I3Context &context) : I3ConditionalModule(context)
    {
        AddParameter("PythiaCommands",
                     "Commands to pass to Pythia.",
                     pythiaCommands_);
        pythiaCommandPath_ = "./pythia-wrapper/resources/cmnd/neutrinoElectroweak.cmnd";
        AddParameter("PythiaCommandPath",
                     "Path to the .cmnd file with pythia commands.",
                     pythiaCommandPath_);
        AddParameter("Statistics Function", "Function to run statistics based off Pythia events.");
    }

    virtual ~I3PythiaWrapperModule();

    virtual void Configure();

    void DAQ(I3FramePtr frame);

    bool ShouldDoProcess(I3FramePtr frame);

    /// @brief Lets you define a custom statistics function callback that gets called after every Pythia event.
    /// @param stats The statistics callback function, taking a Pythia8 Event as an argument.
    void SetStatsFunction(void (*stats)(Pythia8::Event &))
    {
        stats_ = stats;
    }

    void Finish() override;


private:
    I3PythiaWrapperModule();
    I3PythiaWrapperModule(const I3PythiaWrapperModule &);
    I3PythiaWrapperModule &operator=(const I3PythiaWrapperModule &);

    std::vector<std::string> pythiaCommands_;
    std::string pythiaCommandPath_;
    std::map<int, int> pdgHistogram_;
    void (*stats_)(Pythia8::Event &) = nullptr;
    Pythia8::Pythia pythia_;

};

#endif // I3PYTHIAWRAPPERMODULE_HH
