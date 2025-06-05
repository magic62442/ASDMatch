//
// Created by Qiyan LI on 2024/3/4.
//

#include "command.h"

Command::Command(int argc, char **argv) : CommandParser(argc, argv){
    optionsKey[OptionKeyword::QueryGraphPath] = "-q";
    optionsKey[OptionKeyword::DataGraphPath] = "-d";
    optionsKey[OptionKeyword::ResultPath] = "-r";
    optionsKey[OptionKeyword::BaselineType] = "-a";
    optionsKey[OptionKeyword::MemoryBudget] = "-m";
    optionsKey[OptionKeyword::DPStructPath] = "-s";
    optionsKey[OptionKeyword::IntegerToUse] = "-i";
    optionsKey[OptionKeyword::CountPath] = "-c";
    intOptionValue[OptionKeyword::BaselineType] = 0;
    intOptionValue[OptionKeyword::MemoryBudget] = 0;
    processOptions();
}

void Command::processOptions() {
    optionsValue[OptionKeyword::QueryGraphPath] = getCommandOption(optionsKey[OptionKeyword::QueryGraphPath]);
    optionsValue[OptionKeyword::DataGraphPath] = getCommandOption(optionsKey[OptionKeyword::DataGraphPath]);
    optionsValue[OptionKeyword::ResultPath] = getCommandOption(optionsKey[OptionKeyword::ResultPath]);
    intOptionValue[OptionKeyword::BaselineType] = commandOptionType(optionsKey[OptionKeyword::BaselineType]);
    intOptionValue[OptionKeyword::MemoryBudget] = commandOptionType(optionsKey[OptionKeyword::MemoryBudget]);
    intOptionValue[OptionKeyword::IntegerToUse] = commandOptionType(optionsKey[OptionKeyword::IntegerToUse]);
    optionsValue[OptionKeyword::DPStructPath] = getCommandOption(optionsKey[OptionKeyword::DPStructPath]);
    optionsValue[OptionKeyword::CountPath] = getCommandOption(optionsKey[OptionKeyword::CountPath]);
}