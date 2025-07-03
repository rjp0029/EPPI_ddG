/* Created by the PROTEIN PANT(z) Lab at Auburn University.
 *
 * This file contains the definition of the EPPI::BCProps class. This class is a
 * container of the base, primary, and derived features for a specific Binding
 * Complex */

// Make sure the file is only included in a compiled program one time
#ifndef EPPI_BCProps_Guard
#define EPPI_BCProps_Guard 1

// Make sure the file is being included as expected
#ifndef EPPI_Loading_Status
#error BCProps.h has to be included by EPPI.h
#endif

// Define the class
class EPPI::BCProps {

    // The information stored in the class is private
    private:
        // The name of the binding complex
        string m_name;
        // The base features
        vector<double> m_baseFeatures;
        // The primary features
        vector<double> m_primaryFeatures;
        // The derived features
        vector<double> m_derivedFeatures;

    // Private methods to control class behavior
    private:
        // Assign default variables to the class
        void initialize ();
        // Copy content from another instance of the class
        void copy (const BCProps&);
        // Find the index of a base feature
        size_t base_feature_index (const string&);
        // Find the index of a primary feature
        size_t primary_feature_index (const string&);
        // Construct from a set of interactions
        void interaction_constructor(const string& name, vector<HydrogenBond*> hydrogen_bonds,
            vector<SaltBridge*> salt_bridges, vector<Hydrophobic*> hydrophobic_interactions,
            vector<double> ria_values);
        // A constructor function for calculating the values from the files in a
        // directory
        void directory_constructor (const string&, const string&, const bool);
        // Error check that the primary features have appropriate values
        void error_check_primary ();
        // Calculate the primary features from the base features
        void calculate_primary_features ();
        // Calculate a specific derived feature
        void calculate_derived_feature (const Feature&);
        // Calculate all derived features
        void calculate_all_derived_features (const vector<EPPI::Feature>&);
        // these are the functions to parse the 3 files
        map<string, double> parse_per_res_file(string per_res_rosetta_file);
        vector<double> parse_ria_file(string ria_file);
        void parse_eppi_file(string eppi_file, vector<HydrogenBond*> &hydrogen_bonds, 
            vector<SaltBridge*> &salt_bridges, vector<Hydrophobic*> &hydrophobic_interactions);
        // Construct a class instance from "good string" and "row string"
        // outputs
        void good_str_constructor (const vector<string>&);
        void row_str_constructor (const string&, const vector<string>&);

    // The public interface of the class
    public:
        // A default constructor
        BCProps () {initialize();}
        // Copy construction and assignment
        BCProps (const BCProps& other) {initialize(); copy(other);}
        void operator= (const BCProps& other) {copy(other);}
        // Construct from a directory but only the primary features
        BCProps (const string& name, const string& folder, const bool verbose = false) {
            directory_constructor(name, folder, verbose);
            calculate_primary_features ();
        }
        // Construct from a directory and calculate the derived features
        BCProps (const string& name, const string& folder, const vector<Feature>& features, 
                 const bool verbose = false) {
            directory_constructor(name, folder, verbose);
            calculate_primary_features ();
            calculate_all_derived_features (features);}
        // Constructors for loading a BCProps from a "good string" and a "row
        // string"
        BCProps (const vector<string>& lines) {
            good_str_constructor (lines);
            calculate_primary_features ();}
        BCProps (const string& label, const vector<string>& parts) {
            row_str_constructor (label, parts);
            calculate_primary_features ();}
        BCProps (const vector<string>& lines, const vector<Feature>& features) {
            good_str_constructor(lines);
            calculate_primary_features ();
            calculate_all_derived_features (features);}
        BCProps (const string& label, const vector<string>& parts, 
                 const vector<Feature>& features) {
            row_str_constructor(label, parts);
            calculate_primary_features ();
            calculate_all_derived_features (features);}
        // the base features
        vector<double> base_features () const {return m_baseFeatures;}
        // the primary features
        vector<double> primary_features () const {return m_primaryFeatures;}
        // the derived features
        vector<double> derived_features () const {return m_derivedFeatures;}
        // The number of derived features
        size_t size () const {return m_derivedFeatures.size();}
        // Access to a derived feature
        double operator[] (const size_t) const;
        double feature (const size_t i) const {return operator[](i);}
        // naming functions
        void set_name (const string& name) {m_name = name;}
        string name () const {return m_name;}
        // A string that summarizes the base features of the class in a nice
        // manner
        string good_str () const;
        // A string that summarizes the base features of the class in a single
        // row
        string row_str () const;
        // append to csv file
        void append_to_csv (string& filename);

    // End the class definition
};

// Include header files for the class methods
#define EPPI_BCProps_Loading_Status 1
#include "BCProps/initialize.h"
#include "BCProps/copy.h"
#include "BCProps/indexes.h"
#include "BCProps/constructors.h"
#include "BCProps/error_check.h"
#include "BCProps/calculate_features.h"
#include "BCProps/parse_files.h"
#include "BCProps/str_constructors.h"
#include "BCProps/access.h"
#include "BCProps/str.h"
#undef EPPI_BCProps_Loading_Status

// End the header guard from the start of the file
#endif
