//
// Created by anonymous authors on 2024/3/4.
//

#ifndef IN_MEMORY_JOIN_COMMAND_PARSER_H
#define IN_MEMORY_JOIN_COMMAND_PARSER_H


#include <string>
#include <algorithm>
#include <vector>
class CommandParser {
private:
    std::vector<std::string> tokens_;

public:
    CommandParser(int argc, char **argv);
    std::string getCommandOption(const std::string &option) const;
    bool commandOptionExists(const std::string &option) const;
    int commandOptionType(const std::string &option) const;
};


#endif //IN_MEMORY_JOIN_COMMAND_PARSER_H
