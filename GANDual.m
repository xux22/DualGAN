%% GAN
clear all; close all;
rng('default');
profile on;

%% load data
load('mnistAll.mat')
trainX = preprocess(mnist.train_images); 
testX = preprocess(mnist.test_images); 
testY = mnist.test_labels;

% figure;
% idx = randperm(Ntrain); idx = idx(1:64);
% display_mnist(trainX(idx,:));

%% general parameters
settings.epochs = 50;
settings.eta = 0.0002;
settings.alpha = 0.001; % to check
settings.decrease_const = 0.00001;
settings.minibatches = 60000; % not used
settings.batch_size = 1;


%% set generator model parameters
generator_model.n_output = size(trainX,1);
generator_model.n_input = 100; % noisy generator input (previously n_features)
generator_model.n_hidden1 = 256;
generator_model.n_hidden2 = 512;
generator_model.n_hidden3 = 1024;
generator_model.l1 = 0; % not used
generator_model.l2 = 0.1; % not used

%% set discriminator model parameters
discriminator_model.n_output = 1;
discriminator_model.n_input = size(trainX,1);
discriminator_model.n_hidden1 = 1024;
discriminator_model.n_hidden2 = 512;
discriminator_model.n_hidden3 = 256;
discriminator_model.l1 = 0; % not used
discriminator_model.l2 = 0.1; % not used

GAN_fit(settings, generator_model, discriminator_model, trainX);

%% preprocess
function x = preprocess(x)
x = double(x)/255;
x = (x-.5)/.5;
x = reshape(x,28*28,[]);
end
