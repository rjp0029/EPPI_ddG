/* Created by the PROTEIN PANT(z) Lab at Auburn University.
 *
 * This file contains the definition and implementation of the CLASSIFIER
 * protocol. This protocol is used both for training PANTZ classifiers and for
 * using them to evaluate whether predicted complexes are likely to bind
 * experimentally. */

// Make sure the file is being included as expected.
#ifndef Protocol_Loading_Status
#error CLASSIFICATION.h must be included in a compiled program by Protocols.h
#endif

// Define the class
class PROTOCOL::CLASSIFIER {

    // The information stored in the class is private
    private:
        // The name of the calculations
        string m_name;
        // A summary file of the things the calculations have done
        ofstream m_output;
        // The features the Classifier uses
        vector<CLASSIFY::Feature> m_features;
        // The Classifiers that do the classifying
        vector<CLASSIFY::Classifier> m_classifiers;
        // The Complexes used for either training the classifier or that are
        // being evaluated by the classifier
        vector<vector<CLASSIFY::Complex> > m_complexes;
        // The type of calculations that are being done - training or evaluating
        string m_type;
        // A string of information that the classifier generates during
        // construction that should be included in its output file
        string m_message;
        // For training, two numerical parameters are needed
        size_t m_limit;
        double m_P;
        // Header lines from the complex files
        vector<string> m_headers;

    // Private methods of the class that are used during construction of the
    // class
    private:
        // Construction for training calculations
        void training_construction (Interface&);
        // Construction for evaluation calculations
        //void evaluation_construction (Interface&);
        // How to run training calculations
        void training_run ();
        // How to run evaluation calculations
        //void evaluation_run ();

    // The public interface of the class
    public:
        // The destructor
        ~CLASSIFIER () {if (m_output.is_open()) {m_output.close();}}
        // The class constructor load all necessary information and makes sure
        // everything can run.
        CLASSIFIER (Interface&);
        // Run the actual calculations
        void run ();

    // End the class definition
};

// Implement the methods of the Classifier class.

// The constructor for training classifiers
void PROTOCOL::CLASSIFIER::training_construction (Interface& interface) {
    // In addition to the three commands that were already provided, 5
    // additional commands are needed. They are: 1) the file containing the
    // features, 2) the file containing the positive training data, 3) the file
    // containing the negative training data, 4) the minimum number of complexes
    // to include in a classifier, and 5) the percentage range of the positive
    // training data to use for the scoring functions. Thus, 8 total commands
    // are needed.
    size_t N = interface.commands();
    if (N != 8) {
        stringstream c1; c1 << N;
        string error = "Training CLASSIFIER calculations expect exactly 8 "
                       "user-provided instructions, not " + c1.str() + "\n";
        throw PANTZ_error (error);}
    // Set a string to indicate that this is a classifier calculation
    string protocolType = "CLASSIFIER";
    // Get the 4th command
    size_t index = 3;
    string commandType = interface.command_type (index);
    vector<string> * command = interface.command(index);
    // Make sure it is the feature file
    if (commandType != "feature file") {
        string error = "The fourth command for training CLASSIFIER "
                       "calculations must be 'feature file', not '"
                     + commandType + "'\n";
        throw PANTZ_error (error);}
    // Make sure there is only a single entry in the command
    PROTOCOL::validate_command (command, 1, protocolType, commandType);
    // Store the name of the file the features are in. They're loaded after all
    // commands are processed
    string featureFileName = (*command)[0];
    // The 5th command is the file containing the positive training data
    index = 4;
    commandType = interface.command_type (index);
    command = interface.command(index);
    if (commandType != "positive training data file") {
        string error = "The fifth command for training CLASSIFIER calculations "
                       "must be 'positive training data file', not '"
                     + commandType + "'\n";
        throw PANTZ_error (error);}
    PROTOCOL::validate_command (command, 1, protocolType, commandType);
    string goodFileName = (*command)[0];
    // The 6th command is the negative training data
    index = 5;
    commandType = interface.command_type(index);
    command = interface.command(index);
    if (commandType != "negative training data file") {
        string error = "The sixth command for training CLASSIFIER calculations "
                       "must be 'negative training data file', not '"
                     + commandType + "'\n";
        throw PANTZ_error (error);}
    PROTOCOL::validate_command (command, 1, protocolType, commandType);
    string badFileName = (*command)[0];
    // The 7th command is the minimum number of complexes to allow for analysis
    index = 6;
    commandType = interface.command_type(index);
    command = interface.command(index);
    if (commandType != "minimum group size") {
        string error = "The seventh command for training CLASSIFIER "
                       "calculations must be 'minimum group size', not '"
                     + commandType + "'\n";
        throw PANTZ_error (error);}
    PROTOCOL::validate_command (command, 1, protocolType, commandType);
    string number = (*command)[0];
    // Validate that the value is an integer
    if (!Text::is_integer(number)) {
        string error = "The mininum group size must be an integer.\n";
        throw PANTZ_error (error);}
    // Convert it to a number
    stringstream u1; u1 << number; int value; u1 >> value;
    // Make sure the limit is >= 1
    if (value <= 0) {
        string error = "The minimum group size must be a positive integer.\n";
        throw PANTZ_error (error);}
    m_limit = value;
    // The 8th command is the percentage of the good values to use for scoring
    index = 7;
    commandType = interface.command_type(index);
    command = interface.command(index);
    if (commandType != "percentage range") {
        string error = "The eigth command for training CLASSIFIER calculations "
                       "must be 'percentage range', not '" + commandType + "'\n";
        throw PANTZ_error (error);}
    PROTOCOL::validate_command(command, 1, protocolType, commandType);
    number = (*command)[0];
    // Validate that it is a number
    if (!Text::is_number(number)) {
        string error = "The percentage range must be a number.\n";
        throw PANTZ_error (error);}
    stringstream u2; u2 << number; u2 >> m_P;
    if ((m_P <= 0) || (m_P >= 1)) {
        string error = "The percentage range must be greater than 0 and less "
                       "than 1.\n";
        throw PANTZ_error (error);}
    // At this point, all commands have been loaded. Now the data can be loaded
    CLASSIFY::load_features (m_features, featureFileName);
    m_complexes.resize(2);
    CLASSIFY::load_complexes(m_complexes[0], m_features, goodFileName);
    // Get the header line from the good complex file
    ifstream input1; input1.open(goodFileName.c_str());
    string line; getline(input1, line); input1.close();
    Text::strip(line); m_headers.push_back(line);
    CLASSIFY::load_complexes(m_complexes[1], m_features, badFileName);
    // Get the header line from the bad complex file
    input1.open(badFileName.c_str()); getline(input1, line); input1.close();
    Text::strip(line); m_headers.push_back(line);
    // Create the message that will be included at the start of the summary file
    // for these calculations
    stringstream u3; u3 << m_features.size();
    stringstream u4; u4 << m_complexes[0].size();
    stringstream u5; u5 << m_complexes[1].size();
    m_message = m_name + " Classifier Training Calculations\n"
                "Feature File: " + featureFileName + "\n"
                "Features: " + u3.str() + "\n"
                "Positive Training Data File: " + goodFileName + "\n"
                "Positive Training Data Points: " + u4.str() + "\n"
                "Negative Training Data File: " + badFileName + "\n"
                "Negative Training Data Points: " + u5.str() + "\n"
                "Minimum Group Size: " + (*interface.command(6))[0] + "\n"
                "Percentage Range: " + (*interface.command(7))[0] + "\n\n";
    // End the function
}

/*
// Construction for evaluation calculations
void PROTOCOL::CLASSIFIER::evaluation_construction (Interface& interface) {
    // In addition to the three commands already checked, this needs 3 more
    // commands: the name of the classifier, the folder it is located in, and
    // the file of data to evaluate
    size_t N = interface.commands();
    if (N != 6) {
        stringstream c1; c1 << N;
        string error = "Evaluation CLASSIFIER calculations expect exactly 6 "
                       "user-provided instructions, not " + c1.str() + "\n";
        throw PANTZ_error (error);}
    // A string for error checking command contents
    string protocolType = "CLASSIFIER";
    // Get the 4th command
    size_t index = 3;
    string commandType = interface.command_type(index);
    vector<string> * command = interface.command(index);
    // Make sure it is the classifier name
    if (commandType != "classifier name") {
        string error = "The fourth command for evaluation CLASSIFIER "
                       "calculations must be 'classifier name', not '"
                     + commandType + "'\n";
        throw PANTZ_error (error);}
    // Make sure there is only a single entry
    PROTOCOL::validate_command (command, 1, protocolType, commandType);
    // Store the information
    string classifierName = (*command)[0];
    // The 5th command is the folder the classifier is located in
    index = 4;
    commandType = interface.command_type(index);
    command = interface.command(index);
    if (commandType != "classifier folder") {
        string error = "The fifth command for evaluation CLASSIFIER "
                       "calculations must be 'classifier folder', not '"
                     + commandType + "'\n";
        throw PANTZ_error (error);}
    PROTOCOL::validate_command (command, 1, protocolType, commandType);
    string classifierFolder = (*command)[0];
    // Make sure the folder ends with a slash
    if (!Text::endswith(classifierFolder, '/')) {
        classifierFolder += "/";}
    // The 6th command is the name of the file to evaluate
    index = 5;
    commandType = interface.command_type(index);
    command = interface.command(index);
    if (commandType != "evaluation file") {
        string error = "The sixth command for evaluation CLASSIFIER "
                       "calculations must be 'evaluation file', not '"
                     + commandType + "'\n";
        throw PANTZ_error (error);}
    PROTOCOL::validate_command(command, 1, protocolType, commandType);
    string evaluationFile = (*command)[0];
    // Now that everything has been identified, load the information. Make sure
    // to identify the proper name of the file
    string fileName = classifierFolder + "PANTZ_" + classifierName + "_Classifier.txt";
    CLASSIFY::load_classifier(m_classifiers, m_features, fileName);
    // Now load the complexes
    m_complexes.resize(1);
    CLASSIFY::load_complexes(m_complexes[0], m_features, evaluationFile);
    // Create a message that will be included in the output file
    stringstream u1; u1 << m_features.size();
    stringstream u2; u2 << m_classifiers.size();
    stringstream u3; u3 << m_complexes[0].size();
    m_message = m_name + " Evaluation Classifier Calculations\n"
                "Classifier Name: " + classifierName + "\n"
                "Classifier Folder: " + classifierFolder + "\n"
                "File for Evaluation: " + evaluationFile + "\n"
                "Features: " + u1.str() + "\n"
                "Classification Categories: " + u2.str() + "\n"
                "Data Points to Evaluate: " + u3.str() + "\n\n";
    // End the function
}
*/

// Do training calculations
void PROTOCOL::CLASSIFIER::training_run () {
    // The name of the summary file for the calculations
    string summaryFileName = m_name + "_Summary.txt";
    m_output.open(summaryFileName.c_str());
    // Write the summary message to the file
    m_output << m_message;
    // Make a new message listing when the calculations started
    m_message = "Training calculations started on " + METHODS::time_stamp()
              + "\n";
    m_output << m_message << endl;
    // Find the bounds for the features
    /*vector<double> values; values.resize(m_complexes[0].size());
    for(size_t i=0; i<m_features.size(); ++i) {
        for(size_t j=0; j<m_complexes[0].size(); ++j) {
            values[j] = m_complexes[0][j][i];}
        m_features[i].find_bounds(values);}
    // Create an empty vector of Score objects that are the initial
    // qualifications to reach a Classifier
    vector<CLASSIFY::Score> qualifications;
    // Run the recursive training function
    CLASSIFY::recursive_training (m_classifiers, m_complexes[0], 
              m_complexes[1], qualifications, m_features, m_P, m_limit, 0);
    */
    // Make a classifier
    CLASSIFY::Classifier result (m_features, m_complexes[0], m_complexes[1],
                                 m_headers[0], m_headers[1]);
    // Indicate that the classifier training has ended
    m_message = "Training calculations ended on " + METHODS::time_stamp() + "\n";
    m_output << m_message << endl;
    /*
    // Include summaries of the classifiers
    for(size_t i=0; i<m_classifiers.size(); ++i) {
        m_output << m_classifiers[i].str();}
    */
    // Close the output file
    m_output.close();
    // Create the machine-loadable classifier file
    string classifierFileName = "PANTZ_" + m_name + "_Classifier.txt";
    m_output.open(classifierFileName.c_str());
    /*
    // Output the features
    m_output << "FEATURES: " << m_features.size() << "\n";
    for(size_t i=0; i<m_features.size(); ++i) {
        m_output << m_features[i].name() << "\n";}
    // Output the machine versions of the classifiers
    for(size_t i=0; i<m_classifiers.size(); ++i) {
        m_output << m_classifiers[i].machine_str();}
    */
    result.write(m_output);
    // Close the file
    m_output.close();
}

/*
// Do evaluation calculations
void PROTOCOL::CLASSIFIER::evaluation_run () {
    // The name of the file that summarizes the results
    string fileName = m_name + "_Classifications.txt";
    m_output.open(fileName.c_str());
    m_output << m_message;
    // State when calculations happen
    m_message = "Classification calculations started on " + METHODS::time_stamp()
              + "\n";
    m_output << m_message << endl;
    // Loop through the complexes
    for(size_t i=0; i<m_complexes[0].size(); ++i) {
        // Set the message to an empty line
        m_message = "";
        // Classify the complex
        double result = CLASSIFY::evaluate_complex (m_complexes[0][i],
                        m_classifiers, m_message);
        // Output the message
        m_output << m_message;}
    // State when the calculations finished
    m_message = "Classification calculations ended on " + METHODS::time_stamp()
              + "\n";
    m_output << m_message;
    // Close the output file
    m_output.close();
    // End the function
}
*/

// The constructor of the class
PROTOCOL::CLASSIFIER::CLASSIFIER (Interface& interface) {
    // Validate that the interface has at least 3 commands - the type
    // (classification), the name, and the kind of classification
    size_t N = interface.commands();
    if (N < 3) {
        string error = "Algorithm error: CLASSIFIER calculations being "
                       "initialized with insufficient commands.\n";
        throw PANTZ_error (error);}
    // Indicate that this is a CLASSIFIER calculation
    string protocolType = "CLASSIFIER";
    // Get the first command
    size_t index = 0;
    string commandType = interface.command_type(index);
    vector<string> * command = interface.command(index);
    // If the command type is wrong, throw an error
    if (commandType != "calculation type") {
        string error = "Algorithm error: the first command for CLASSIFIER "
                       "calculations must be 'calculation type', not '"
                     + commandType + "'\n";
        throw PANTZ_error (error);}
    // Validate that there is one entry
    PROTOCOL::validate_command (command, 1, protocolType, commandType);
    // Convert it to lower case and confirm that it is 'classifier'
    string text = (*command)[0]; Text::lower(text);
    if (text != "classifier") {
        string error = "Algorithm error: CLASSIFIER object being constructed "
                       "for " + (*command)[0] + " calculations.\n";
        throw PANTZ_error (error);}
    // The second command must be the name of the calculations
    index = 1;
    commandType = interface.command_type(index);
    command = interface.command(index);
    if (commandType != "calculation name") {
        string error = "The second command in CLASSIFIER calculations must be "
                       "'calculation name', not " + commandType + "\n";
        throw PANTZ_error(error);}
    PROTOCOL::validate_command (command, 1, protocolType, commandType);
    m_name = (*command)[0];
    // The third command must by the kind of classification calculation
    index = 2;
    commandType = interface.command_type(index);
    command = interface.command(index);
    if (commandType != "kind of calculation"){
        string error = "The third command in CLASSIFIER calculations must be "
                       "'kind of calculation', not " + commandType + "\n";
        throw PANTZ_error (error);}
    PROTOCOL::validate_command (command, 1, protocolType, commandType);
    m_type = (*command)[0];
    // Convert the type to lower case
    Text::lower(m_type);
    // Finish construction using the specialized constructors
    if (m_type == "training") {training_construction(interface);}
    //else if (m_type == "evaluation") {evaluation_construction(interface);}
    // Otherwise, throw an error
    else {
        string error = "'" + m_type + "' is not a recognized kind of "
                       "calculation for CLASSIFIER calculations.\n";
        throw PANTZ_error (error);}
    // End the constructor
}

// Run the classifier calculations
void PROTOCOL::CLASSIFIER::run () {
    // Do this based on the type of calculation
    if (m_type == "training") {training_run();}
    //else if (m_type == "evaluation") {evaluation_run();}
    // Otherwise, throw an error
    else {
        string error = "Algorithm error: the run function of the CLASSIFIER "
                       "class does not support " + m_type + " calculations\n";
        throw PANTZ_error (error);}
}
