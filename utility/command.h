//
// Created by anonymous authors on 2024/3/4.
//

#ifndef IN_MEMORY_JOIN_COMMAND_H
#define IN_MEMORY_JOIN_COMMAND_H

#include "command_parser.h"
#include <iostream>
#include <map>

enum OptionKeyword {
    QueryGraphPath = 1,     // -q, the query graph file path
    DataGraphPath = 2,      // -d, the data graph file path
    ResultPath = 3,         // -r, the result file path
    BaselineType = 4,        // -a, the baseline type
    MemoryBudget = 5,         // -m, the memory budget
    DPStructPath = 6,         // -s, the dp structure file path
    IntegerToUse = 7,         // -i, some integer to use
};

class Command : public CommandParser {
private:
    std::map<OptionKeyword, std::string> optionsKey;
    std::map<OptionKeyword, std::string> optionsValue;
    std::map<OptionKeyword, bool> booleanOptionValue;
    std::map<OptionKeyword, int> intOptionValue;

private:
    void processOptions();

public:
    Command(int argc, char **argv);

    std::string getQueryGraphPath() {
        return optionsValue[OptionKeyword::QueryGraphPath];
    }

    std::string getDataGraphPath() {
        return optionsValue[OptionKeyword::DataGraphPath];
    }

    std::string getResultPath() {
        return optionsValue[OptionKeyword::ResultPath];
    }

    int getBaselineType() {
        return intOptionValue[OptionKeyword::BaselineType];
    }

    int getMemoryBudget() {
        return intOptionValue[OptionKeyword::MemoryBudget];
    }

    int getInteger() {
        return intOptionValue[OptionKeyword::IntegerToUse];
    }

    std::string getDPStructPath() {
        return optionsValue[OptionKeyword::DPStructPath];
    }
};


#endif //IN_MEMORY_JOIN_COMMAND_H
