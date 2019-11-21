// -*- c++ -*-
#ifndef EventSelection_KerasModel_h
#define EventSelection_KerasModel_h

#include <iostream>
#include <fstream>

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
  virtual std::vector<float> predict(std::vector<float> test_input){
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
    float f_weight;
    // Read the output weights for each input nodes
    char tmp_char;
    for (unsigned int i = 0; i < input_nodes; i++) {
      //std::cout<<i<<" / "<<input_nodes<<std::endl;
      fin >> tmp_char; // '['
      std::vector<float> v_weights;
      for (unsigned int j = 0; j < output_nodes; j++) {
	// Read values and store in vector
	fin >> f_weight;
	v_weights.push_back(f_weight);
      }
      fin >> tmp_char; // ']'
      layer_weights.push_back(v_weights); //layer_weights: vector with dimensions input_nodes x output_nodes
    }
    // Read and store bias values
    fin >> tmp_char; // '['                                                                                                                                                                  
    for (unsigned int k = 0; k < output_nodes; k++) {
      fin >> f_weight;
      bias.push_back(f_weight);
    }
    fin >> tmp_char; // ']' 
  }

  //=== LayerDense: Compute the output of a Dense layer
  std::vector<float> predict(std::vector<float> test_input)
    {
      // Returns a vector with size output_nodes.
      std::vector<float> out(output_nodes);
      float weighted_term = 0;
      for (size_t i = 0; i < output_nodes; i++) {
	weighted_term = 0;
	for (size_t j = 0; j < input_nodes; j++) {
	  //weighted_term_i = Sum_j (input_j * W_ji)
	  weighted_term += (test_input[j] * layer_weights[j][i]);
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
  std::vector<float> predict(std::vector<float> test_input)
    {
      //List of activation functions: https://keras.io/activations/
      if(activation_type == "softmax") {
	// y_i = exp(x_i)/(Sum_j exp(x_j))
	float sum = 0.0;
        for(unsigned int k = 0; k < test_input.size(); ++k) {
	  test_input[k] = exp(test_input[k]);
	  sum += test_input[k];
        }
        for(unsigned int k = 0; k < test_input.size(); ++k) {
	  test_input[k] /= sum;
        }
      }
      else if(activation_type == "softplus") {
	// y_i = log(1+exp(x_i))
	for (unsigned int k = 0; k < test_input.size(); ++k) {	                                                                                                        
	  test_input[k] = log1p(exp(test_input[k]));  // log1p = (ln(1+x))
	}
      }
      else if(activation_type == "softsign") {
	// y_i = x_i/(1 + abs(x_i))
	for (unsigned int k = 0; k < test_input.size(); ++k) {
	  test_input[k] = test_input[k]/(1+abs(test_input[k]));
	}
      }
      else if(activation_type == "relu") {
	// y_i = max(0,x_i) (With default values!)
	for (unsigned int i = 0; i < test_input.size(); i++) {
	  if(test_input[i] < 0) {
	    test_input[i] = 0;
	  }
	}
      }
      else if (activation_type == "tanh") {
	// y_i = tanh(x_i) = (exp(x_i) - exp(-x_i)) / (exp(x_i) + exp(-x_i))
	for(unsigned int k = 0; k < test_input.size(); ++k) {
	  float denominator = exp(test_input[k]) + exp(-test_input[k]);
	  float numerator   = exp(test_input[k]) - exp(-test_input[k]);
	  test_input[k] = numerator / denominator;
        }
      }
      else if (activation_type == "sigmoid") {
	// y_i = 1/(1 + exp(-x_i))
	for(unsigned int k = 0; k < test_input.size(); ++k) {
	  float denominator = 1 + exp(-(test_input[k]));
	  test_input[k] = 1/denominator;
        }
      }
      else if (activation_type == "hard_sigmoid") {
	// -> y_i = 0 (if x_i < -2.5)
	// -> y_i = 1 (if x_i > +2.5)
	// -> y_i = 0.2*x_i + 0.5 (if -2.5 <= x_i <= 2.5)
	for(unsigned int k = 0; k < test_input.size(); ++k) {
	  
	  if (test_input[k] < -2.5){
	    test_input[k] = 0;
	  }
	  else if (test_input[k] > 2.5){
	    test_input[k] = 1;
	  }
	  else{
	    test_input[k] = 0.2*test_input[k] + 0.5;
	  }
	}
      }      
      else if (activation_type == "linear") {
	// y_i = x_i
	return test_input;
      }
      else {
	std::cout << "Activation " << activation_type << " not defined!" << std::endl;
	return test_input;
      }
      /*
      else if (activation_type == "elu") {
	// -> y_i = x_i (if x_i > 0)
	// -> y_i = alpha * (exp(x) -1) (if x_i < 0)
	float alpha = 1; //Default: should be fixed if elu is used!
	for (unsigned int i = 0; i < test_input.size(); i++) {
	  if(test_input[i] < 0) {
	    test_input[i] = alpha*(exp(test_input[i]-1));
	  }
	}
      }
      else if (activation_type == "selu") {
	// y_i = scale*elu(x, alpha)
	float alpha = 1; //Default: should be fixed!
	float scale = 1; //Is that given as parameter? To be fixed! 
	for (unsigned int i = 0; i < test_input.size(); i++) {
	  if(test_input[i] < 0) {
	    test_input[i] = scale*alpha*(exp(test_input[i]-1));
	  }
	  else{
	    test_input[i] = scale*test_input[i];
	  }			
	}
      }
      */
      return test_input;
    }
};

//================================================================================================
// Keras Model class
//================================================================================================
class KerasModel {
 public:
  unsigned int input_nodes();
  unsigned int output_nodes();
	
  KerasModel(){}
  ~KerasModel()
    {
      for (unsigned int i = 0; i < layers.size(); i++) {
	delete layers[i];       // deallocate memory                                                                                                                                                
      }
    }
	
  //=== Keras Model: predict output
  std::vector<float> predict(std::vector<float> test_input){
    std::vector<float> response;
    // Iterate over each layer
    for (unsigned int i = 0; i < n_layers; i++) {
      if (0) std::cout << "Processing layer to compute output " << layers[i]->layer_name << std::endl;
      // Compute output (based on the type of the layer)
      response = layers[i]->predict(test_input);
      // The output of each layer is the input of the nextt layer
      test_input = response;
    }
    return response;
  }
  
  //Number of layers
  unsigned int n_layers, n_inputs;
  std::vector<std::string> inputs;
  //Vector of layers
  std::vector<Layer *> layers;
  
  //=== Keras Model: Load weights of each layer
  void load_weights(std::string &input_fname)
  {
    if (0) std::cout << "=== KerasModel: Reading weights from file " << input_fname << std::endl;
    std::ifstream fin(input_fname, std::ifstream::in);
    std::string layer_type, input_name;
    std::string  _str; //temporary variable
    
    if(fin.is_open()) {
      // Get list of inputs
      fin >> _str >> n_inputs;
      for (unsigned int input_index = 0; input_index < n_inputs; ++input_index) {
	fin >> input_name;
	inputs.push_back(input_name);
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
	// if none of above case is true, means layer not-defined
	if(layer == 0L) {
	  std::cout << "=== KerasModel: Layer of type " << layer_type << " is not defined" << std::endl;
	  return;
	}
      
	// Load model weights
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
