//
//  SimOptions.cpp
//
//  Christopher Batty, Fang Da 2014
//
//

#include <assert.h>
#include <fstream>
#include <sstream>
#include "SimOptions.h"

std::map<std::string, Options::Option> Options::s_options;

void Options::addStringOption(const std::string & key, const std::string & default_value)
{
    assert(s_options.find(key) == s_options.end()); // verify this option doesn't already exit
    
    Option o;
    o.key = key;
    o.type = STRING;
    o.str_value = default_value;
    
    s_options[key] = o;
}

void Options::addIntegerOption(const std::string & key, int default_value)
{
    assert(s_options.find(key) == s_options.end()); // verify this option doesn't already exit
    
    Option o;
    o.key = key;
    o.type = INTEGER;
    o.int_value = default_value;
    
    s_options[key] = o;
}

void Options::addDoubleOption(const std::string & key, double default_value)
{
    assert(s_options.find(key) == s_options.end()); // verify this option doesn't already exit
    
    Option o;
    o.key = key;
    o.type = DOUBLE;
    o.double_value = default_value;
    
    s_options[key] = o;
}

void Options::addBooleanOption(const std::string & key, bool default_value)
{
    assert(s_options.find(key) == s_options.end()); // verify this option doesn't already exit
    
    Option o;
    o.key = key;
    o.type = BOOLEAN;
    o.bool_value = default_value;
    
    s_options[key] = o;
}


bool Options::parseOptionFile(const std::string & file, bool verbose)
{
    std::ifstream fin(file.c_str());
    
    std::string line;
    while (!fin.eof())
    {
        std::getline(fin, line);
        std::stringstream ss(line);
        
        std::string key;
        ss >> key;
        if (key == "#" || key == "" || ss.eof())    // skip comment lines and empty lines
            continue;

        std::map<std::string, Option>::iterator i = s_options.find(key);
        assert(i != s_options.end());
        
        switch (i->second.type)
        {
            case STRING:
                ss >> i->second.str_value;
                break;
            case INTEGER:
                ss >> i->second.int_value;
                break;
            case DOUBLE:
                ss >> i->second.double_value;
                break;
            case BOOLEAN:
                ss >> i->second.bool_value;
                break;
            default:
                assert(!"Unexpected option type");
                break;
        }
    }

    if (verbose)
    {
        for (std::map<std::string, Option>::iterator i = s_options.begin(); i != s_options.end(); i++)
        {
            std::cout << "option " << i->first << " = ";
            switch (i->second.type)
            {
                case STRING:
                    std::cout << i->second.str_value;
                    break;
                case INTEGER:
                    std::cout << i->second.int_value;
                    break;
                case DOUBLE:
                    std::cout << i->second.double_value;
                    break;
                case BOOLEAN:
                    std::cout << i->second.bool_value;
                    break;
                default:
                    assert(!"Unexpected option type");
                    break;
            }
            std::cout << std::endl;
        }
    }
    
    return false;
}

const std::string & Options::strValue(const std::string & key)
{
    assert(s_options.find(key) != s_options.end()); // verify this option exists
    assert(s_options[key].type == STRING);          // verify this option has the correct type
    return s_options[key].str_value;
}

int Options::intValue(const std::string & key)
{
    assert(s_options.find(key) != s_options.end()); // verify this option exists
    assert(s_options[key].type == INTEGER);         // verify this option has the correct type
    return s_options[key].int_value;
}

double Options::doubleValue(const std::string & key)
{
    assert(s_options.find(key) != s_options.end()); // verify this option exists
    assert(s_options[key].type == DOUBLE);          // verify this option has the correct type
    return s_options[key].double_value;
}

bool Options::boolValue(const std::string & key)
{
    assert(s_options.find(key) != s_options.end()); // verify this option exists
    assert(s_options[key].type == BOOLEAN);         // verify this option has the correct type
    return s_options[key].bool_value;
}

