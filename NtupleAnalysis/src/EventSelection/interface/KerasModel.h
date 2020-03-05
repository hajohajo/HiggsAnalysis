// -*- c++ -*-
#ifndef EventSelection_KerasModel_h
#define EventSelection_KerasModel_h

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

class KerasModel;
class Layer;
class LayerDense;
class LayerActivation;

//================================================================================================
// Layer class
//================================================================================================
class Layer {
 public:
  unsigned int layer_id;
  std::string layer_name;

  //Constructor: Sets layer name (Dense / Activation)
  Layer(std::string name) : layer_name(name) {};
  virtual ~Layer() { };
  
  virtual void load_weights(std::ifstream &input_fname) { };
  //Computes the output value of the NN model
  virtual std::vector<float> predict(std::vector<float> x_inputs){
    std::vector<float> response{-1};
    return response; //Default value = -1 (Output should range from 0 to 1)
  }
  //Returns layer name
  std::string getLayerName() { return layer_name; }
  };

//================================================================================================
// Dense Layer class
//================================================================================================
class LayerDense: public Layer {
private:
  void GetWeightsFromArray(std::ifstream &fin, std::vector<float> &v, unsigned int n_weights){	
    //tmp variables
    std::string _str;
    float f_weight;    
    fin >> _str; // '['
    for (unsigned int j = 0; j < n_weights; j++) {
      fin >> f_weight;
      v.push_back(f_weight);
    }
    fin >> _str; // ']'
  }       
  
public:
  // Number of input and output nodes
  unsigned int input_nodes, output_nodes;
  std::vector<std::vector<float> > layer_weights;
  std::vector<float> bias;
  LayerDense() : Layer("Dense") {};

  //=== LayerDense: Load weights
  void load_weights(std::ifstream &fin)
  {
    if (0) std::cout<<"=== LayerDense: Load weights"<<std::endl;
    fin >> input_nodes >> output_nodes;
    if (0) std::cout<<"=== LayerDense: Load weights - input, output"<<input_nodes<<" "<<output_nodes<<std::endl;
    // Read the output weights for each input nodes
    for (unsigned int i = 0; i < input_nodes; i++) {
      std::vector<float> v_weights;
      GetWeightsFromArray(fin, v_weights, output_nodes);
      // layer_weights: vector with dimensions input_nodes x output_nodes
      layer_weights.push_back(v_weights); 
    }
    // Read and store bias values
    GetWeightsFromArray(fin, bias, output_nodes);
  }

  //=== LayerDense: Compute the output of a Dense layer
  std::vector<float> predict(std::vector<float> x_inputs)
    {
      // Returns a vector with size output_nodes.
      std::vector<float> out(output_nodes);
      float weighted_term = 0;
      for (size_t i = 0; i < output_nodes; i++) {
	weighted_term = 0;
	for (size_t j = 0; j < input_nodes; j++) {
	  //weighted_term_i = Sum_j (input_j * W_ji)
	  weighted_term += (x_inputs[j] * layer_weights[j][i]);
	}
	//Output_i = weighted_term_i + bias_i (shifts the activation function)
	out[i] = weighted_term + bias[i];
      }
      return out;
    }
};

//================================================================================================
// Activation Layer class
//================================================================================================
class LayerActivation: public Layer {
 public:
  std::string activation_type;

  LayerActivation() : Layer("Activation") { };

  //=== LayerActivation: Load weights
  void load_weights(std::ifstream &fin)
  {
    fin >> activation_type;
    if (0) std::cout<<"LayerActivation: Activation type "<<activation_type<<std::endl;
  }

  //=== LayerActivation: Compute Output
  std::vector<float> predict(std::vector<float> x_inputs)
    {
      //List of activation functions: https://keras.io/activations/
      if(activation_type == "softmax") {
	// y_i = exp(x_i)/(Sum_j exp(x_j))
	float sum = 0.0;
        for(unsigned int k = 0; k < x_inputs.size(); ++k) {
	  x_inputs[k] = exp(x_inputs[k]);
	  sum += x_inputs[k];
        }
        for(unsigned int k = 0; k < x_inputs.size(); ++k) {
	  x_inputs[k] /= sum;
        }
      }
      else if(activation_type == "softplus") {
	// y_i = log(1+exp(x_i))
	for (unsigned int k = 0; k < x_inputs.size(); ++k) {	                                                                                                        
	  x_inputs[k] = log1p(exp(x_inputs[k]));  // log1p = (ln(1+x))
	}
      }
      else if(activation_type == "softsign") {
	// y_i = x_i/(1 + abs(x_i))
	for (unsigned int k = 0; k < x_inputs.size(); ++k) {
	  x_inputs[k] = x_inputs[k]/(1+abs(x_inputs[k]));
	}
      }
      else if(activation_type == "relu") {
	// y_i = max(0,x_i) (With default values!)
	for (unsigned int i = 0; i < x_inputs.size(); i++) {
	  if(x_inputs[i] < 0) {
	    x_inputs[i] = 0;
	  }
	}
      }
      else if (activation_type == "tanh") {
	// y_i = tanh(x_i) = (exp(x_i) - exp(-x_i)) / (exp(x_i) + exp(-x_i))
	for(unsigned int k = 0; k < x_inputs.size(); ++k) {
	  float denominator = exp(x_inputs[k]) + exp(-x_inputs[k]);
	  float numerator   = exp(x_inputs[k]) - exp(-x_inputs[k]);
	  x_inputs[k] = numerator / denominator;
        }
      }
      else if (activation_type == "sigmoid") {
	// y_i = 1/(1 + exp(-x_i))
	for(unsigned int k = 0; k < x_inputs.size(); ++k) {
	  float denominator = 1 + exp(-(x_inputs[k]));
	  x_inputs[k] = 1/denominator;
        }
      }
      else if (activation_type == "hard_sigmoid") {
	// -> y_i = 0 (if x_i < -2.5)
	// -> y_i = 1 (if x_i > +2.5)
	// -> y_i = 0.2*x_i + 0.5 (if -2.5 <= x_i <= 2.5)
	for(unsigned int k = 0; k < x_inputs.size(); ++k) {
	  
	  if (x_inputs[k] < -2.5){
	    x_inputs[k] = 0;
	  }
	  else if (x_inputs[k] > 2.5){
	    x_inputs[k] = 1;
	  }
	  else{
	    x_inputs[k] = 0.2*x_inputs[k] + 0.5;
	  }
	}
      }      
      else if (activation_type == "linear") {
	// y_i = x_i
	return x_inputs;
      }
      else {
	std::cout << "Activation " << activation_type << " not defined!" << std::endl;
	return x_inputs;
      }
      /*
      else if (activation_type == "elu") {
	// -> y_i = x_i (if x_i > 0)
	// -> y_i = alpha * (exp(x) -1) (if x_i < 0)
	float alpha = 1; //Default: should be fixed if elu is used!
	for (unsigned int i = 0; i < x_inputs.size(); i++) {
	  if(x_inputs[i] < 0) {
	    x_inputs[i] = alpha*(exp(x_inputs[i]-1));
	  }
	}
      }
      else if (activation_type == "selu") {
	// y_i = scale*elu(x, alpha)
	float alpha = 1; //Default: should be fixed!
	float scale = 1; //Is that given as parameter? To be fixed! 
	for (unsigned int i = 0; i < x_inputs.size(); i++) {
	  if(x_inputs[i] < 0) {
	    x_inputs[i] = scale*alpha*(exp(x_inputs[i]-1));
	  }
	  else{
	    x_inputs[i] = scale*x_inputs[i];
	  }			
	}
      }
      */
      return x_inputs;
    }
};

//================================================================================================
// Keras Model class
//================================================================================================
class KerasModel {
private:  
  void GetWeightsFromArray(std::ifstream &fin, std::vector<float> &v, unsigned int n_inputs){	
    //tmp variables
    std::string _str;
    float f_weight;    
    
    fin >> _str; // '['
    for (unsigned int j = 0; j < n_inputs; j++) {
      fin >> f_weight;
      v.push_back(f_weight);
    }
    fin >> _str; // ']'
  }       

  void ReadAndStoreScalerAttributes(std::ifstream &fin, std::vector<float> &mu, std::vector<float> &scale, bool &isSparse, unsigned int n_inputs){
    std::string _str, type;
    for (int i = 0; i < 2; i++){
      fin >> type; //Attribute Type (mu, scale)
      if (0) std::cout<<"=== KerasModel: Attribute type "<< type<<std::endl;
      if (type == "scale"){
	GetWeightsFromArray(fin, scale, n_inputs);
      }
      else{
	// mu = mean (Standard), min (MinMax), center (Robust)
	GetWeightsFromArray(fin, mu, n_inputs);
      }
    }
    // isSparse?    
    fin >> _str >> isSparse;
    if (0) std::cout<<"isSparse? "<<isSparse<<std::endl;
  }  

  void TransformVariables(std::vector<float> &x_inputs, std::string scaler_type,  std::vector<float> mu,  std::vector<float> scale, bool isSparse){
    // Preprocessing: https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/preprocessing/_data.py
    
    //Sanity check
    if (x_inputs.size() != mu.size() || x_inputs.size() != scale.size()){
      std::cout<<"=== KerasModel:TransformVariables: Scaler attribute vectors are not of the same size as the input variables!"<<std::endl;
      return;
    }

    if (scaler_type == "Standard"){
      // X_transformed = (X - mean)/scale
      for (unsigned int i = 0; i < x_inputs.size(); i++){
	x_inputs.at(i) = (x_inputs.at(i) - mu.at(i))/scale.at(i);
      }	     
    }
    else if (scaler_type == "MinMax"){
      //X_transformed = X*scale_ + min_
      for (unsigned int i = 0; i < x_inputs.size(); i++){
	x_inputs.at(i) = x_inputs.at(i)*scale.at(i) + mu.at(i);
      }
    }
    else if (scaler_type == "Robust"){
      //X_transformed = (X-center_)/scale_
      for (unsigned int i = 0; i < x_inputs.size(); i++){
	x_inputs.at(i) = (x_inputs.at(i) - mu.at(i))/scale.at(i);
      }
    }
    else{
      std::cout << "Scaler type " << scaler_type << " not defined!" << std::endl;
      return;
    }
    return;
  }
  
public:
  unsigned int input_nodes();
  unsigned int output_nodes();

  // Number of layers
  unsigned int n_layers, n_inputs;
  std::vector<std::string> inputs;
  // Vector of layers
  std::vector<Layer *> layers;

  // Scaler (input preprocessing)
  std::string scaler_type = "None";
  std::vector<float> mu, scale;
  bool isSparse = false;
	
  KerasModel(){}
  ~KerasModel()
    {
      for (unsigned int i = 0; i < layers.size(); i++) {
	delete layers[i]; // deallocate memory                                                                                                                                                
      }
    }

  std::vector<float> predict(std::vector<float> x_inputs){
    std::vector<float> response;
    // Variable preprocessing  

    if (0){    
      std::cout<<"Inputs before transformation"<<std::endl;
      for (unsigned int i = 0; i < x_inputs.size(); i++) {
	std::cout<<i<<" "<<x_inputs.at(i)<<std::endl;
      }
    }

    if (scaler_type != "None"){
      //Standardization
      TransformVariables(x_inputs, scaler_type, mu, scale, isSparse);
    }

    if (0){
      std::cout<<"Inputs after transformation"<<std::endl;
      for (unsigned int i = 0; i < x_inputs.size(); i++) {
	std::cout<<i<<" "<<x_inputs.at(i)<<std::endl;
      }
    }
    
    // Iterate over each layer
    for (unsigned int i = 0; i < n_layers; i++) {
      //if (0) std::cout << "Processing layer to compute output " << layers[i]->layer_name << std::endl;
      // Compute output (based on the layer type)
      response = layers[i]->predict(x_inputs);
      // The output of each layer is the input of the next layer
      x_inputs = response;
    }
    return response;
  }
  
  //=== Keras Model: Load weights of each layer
  void load_weights(std::string &input_fname)
  {
    if (0) std::cout << "=== KerasModel: Reading weights from file " << input_fname << std::endl;
    std::ifstream fin(input_fname, std::ifstream::in);
    std::string layer_type, input_name;
    std::string  _str; //temporary variable

    if(fin.is_open()) {

      // Get number of inputs
      fin >> _str >> n_inputs;
      // Get list of inputs
      for (unsigned int input_index = 0; input_index < n_inputs; ++input_index) {
	fin >> input_name;
	inputs.push_back(input_name);
      }

      //Get scaler type
      fin >> _str >> scaler_type;

      if (scaler_type != "None"){      
	// Read and store the scaler attributes
	ReadAndStoreScalerAttributes(fin, mu, scale, isSparse, n_inputs);
	
	//Check if training datast is sparse
	if (isSparse){
	  std::cout << "Training dataset isSparse! Features cannot be transformed!"<<std::endl;
	  return;
	}
      }

      // Get number of layers
      fin >> _str >> n_layers;
      if (0) std::cout<<"=== KerasModel: number of layers "<<n_layers<<std::endl;
      
      // Iterate over each layer
      for (unsigned int layer_index = 0; layer_index < n_layers; ++layer_index) {
	// Get layer type
	fin >> layer_type;
	if (0) std::cout << "=== KerasModel: iterating over each layer "<<layer_index << " "<<layer_type<<std::endl;
	
	// Pointer to layer
	Layer *layer = 0L;	

	if (layer_type == "Dense") {
	  layer = new LayerDense();
	}
	else if(layer_type == "Activation") {
	  layer = new LayerActivation();
	}
	// if none of above cases is true, means layer not-defined
	if(layer == 0L) {
	  std::cout << "=== KerasModel: Layer of type " << layer_type << " is not defined" << std::endl;
	  return;
	}
      
	// Load model weights (based on the layer type)
	layer->load_weights(fin);
	// Store layer
	layers.push_back(layer);
      }
    } //if(fin.is_open()) {
    fin.close();	  
  }// void load_weights(std::string &input_fname)...{

  std::vector<std::string> GetInputNames(){
    return inputs;
  }
  
};

#endif
