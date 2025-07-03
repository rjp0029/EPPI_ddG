/* Created by the PROTEIN PANT(z) Lab at Auburn University.
 *
 * This file contains the definition of the EPPI::Feature class. This class is
 * primarily labelling information: the name of the feature, its index in a list
 * of features being used in an analysis, and information on how to calculate
 * the feature from the Primary Features of the EPPI namespace. */

// Make sure the file is only included in a compiled program one time
#ifndef EPPI_Feature_Guard
#define EPPI_Feature_Guard 1

// Make sure the file is being included as expected
#ifndef EPPI_Loading_Status
#error Feature.h has to be included by EPPI.h
#endif

// Define the class
class EPPI::Feature {

    // The information stored in the class is private
    private:
        // The name of the feature
        string m_name;
        // It's index in a list of features
        size_t m_index;
        // And information about the primary features it is calculated from
        vector<string> m_primary;
        vector<char> m_operations;

    // Private methods to control class behavior
    private:
        // Assign default variables to the class
        void initialize ();
        // Copy content from another instance of the class
        void copy (const Feature&);

    // The public interface of the class
    public:
        // A default constructor
        Feature () {initialize();}
        // Copy construction and assignment
        Feature (const Feature& other) {initialize(); copy(other);}
        void operator= (const Feature& other) {copy(other);}
        // Construct a feature
        Feature (const string&, const size_t);
        // Access to feature information
        string name () const {return m_name;}
        size_t index () const {return m_index;}
        size_t primary_features () const {return m_primary.size();}
        string primary_feature (const size_t) const;
        char operation (const size_t) const;

    // End the Feature class definition
};

// Load the feature header files
#define EPPI_Feature_Loading_Status 1
#include "Feature/initialize.h"
#include "Feature/copy.h"
#include "Feature/constructor.h"
#include "Feature/access.h"

// End the header guard from the start of the file
#endif
