


import sympy
import sklearn.datasets
import sklearn.feature_extraction.text
import umap
import umap.plot
import matplotlib.pyplot as plt
import sparse


import pandas as pd
import numpy as np
from scipy.io import mmread
from scipy import sparse
import scipy.sparse
from tensorflow import keras
from tensorflow.keras import layers

from tensorflow.keras.regularizers import l2


import shap


model = keras.Sequential([  
    keras.layers.BatchNormalization(),
    keras.layers.Dense(128),  
    keras.layers.Dropout(0.2),
    keras.layers.BatchNormalization(),
    keras.layers.Dense(32),  
    keras.layers.Dropout(0.2),
    keras.layers.BatchNormalization(),
    keras.layers.Dense(8),  
    keras.layers.Dropout(0.2),
    keras.layers.Dense(1, activation='sigmoid')
])

model = keras.Sequential([  
    keras.layers.BatchNormalization(),
    keras.layers.Dense(128),  
    keras.layers.Dropout(0.2),
    keras.layers.Dense(1, activation='sigmoid')
])

# Input: a single integer per sample
input_layer = keras.layers.Input(shape=(544977,))

# Embedding: maps each integer (0–999) to a 32-dimensional vector
embedding_layer = keras.layers.Embedding(input_dim=544977, output_dim=100)(input_layer)

# Flatten the embedding to feed into dense layers
flattened = keras.layers.Flatten()(embedding_layer)

# Continue with normal dense layers
dense = keras.layers.Dense(64, activation='relu')(flattened)
output = keras.layers.Dense(1, activation='sigmoid')(dense)

model = keras.Model(inputs=input_layer, outputs=output)
model.compile(optimizer='adam', loss='binary_crossentropy', metrics=metrics)



model = keras.Sequential([
  keras.layers.BatchNormalization(),
  keras.layers.Dense(64, activation=keras.activations.get('linear')),
  keras.layers.LeakyReLU(alpha=0.01),
  keras.layers.Dropout(0.5),
  keras.layers.Dense(1, activation='sigmoid')
])


model = keras.Sequential([
    keras.layers.BatchNormalization(),
    keras.layers.Dense(64, activation=keras.activations.get('linear')),
    keras.layers.Dense(64, activation='relu'),
    #keras.layers.LeakyReLU(alpha=0.01),
    keras.layers.Dense(1, activation='sigmoid')
])

model = keras.Sequential([
    keras.Input(shape=(x_train.shape[1],)),
    
    # Optional normalization
    layers.BatchNormalization(),
    
    # First dense block
    layers.Dense(512, activation=None, kernel_regularizer=regularizers.l2(1e-4)),
    layers.LeakyReLU(alpha=0.01),
    layers.Dropout(0.5),

    # Second dense block
    layers.Dense(128, activation=None, kernel_regularizer=regularizers.l2(1e-4)),
    layers.LeakyReLU(alpha=0.01),
    layers.Dropout(0.3),

    # Output layer for binary classification
    layers.Dense(1, activation='sigmoid')
])

callbacks = [
    keras.callbacks.EarlyStopping(
    monitor="val_auc",
    min_delta=0.001,
    patience=3,
    verbose=0,
    mode="auto",
    baseline=None,
    restore_best_weights=True,
    start_from_epoch=20),
    keras.callbacks.ReduceLROnPlateau(patience=3, factor=0.5)
]

metrics = [
    keras.metrics.FalseNegatives(name="fn"),
    keras.metrics.FalsePositives(name="fp"),
    keras.metrics.TrueNegatives(name="tn"),
    keras.metrics.TruePositives(name="tp"),
    keras.metrics.Precision(name="precision"),
    keras.metrics.Recall(name="recall"),
    keras.metrics.AUC(name = "auc"),
    keras.metrics.SensitivityAtSpecificity(0.9, name = "sensitivity_at_spec0.9")
    ]

#class_weight = {0: weight_for_0, 1: weight_for_1}


model.compile(optimizer=keras.optimizers.Adam(learning_rate=1e-3),
              loss='binary_crossentropy',
              metrics=metrics)


model.fit(x_train, 
          y_train, 
          epochs=100,
          callbacks=callbacks,
          validation_data=(x_val, y_val))


model.evaluate(x_test, y_test, verbose=2)

# 6. Make predictions (optional)
predictions = model.predict(x_test[:5]) # Predict on the first 5 test images
predicted_labels = np.argmax(predictions, axis=1)

prediction_scores = model.predict(np.array(dat.todense()))


col_index = res.columns.get_loc('code') + 1
res.insert(col_index, 'contact_score', prediction_scores)

res.to_csv("/oak/stanford/groups/ebutcher/kevin/bm_contact_score_example3.csv")




model.save("/oak/stanford/groups/ebutcher/kevin/model_july13_auc0.82.keras")


model = keras.models.load_model("/oak/stanford/groups/ebutcher/kevin/model_july9_auc0.84_2.keras")

res = pd.read_csv("/oak/stanford/groups/ebutcher/kevin/ds_cons.csv")

explainer = shap.DeepExplainer(model, x_train)

sample_indices = np.random.choice(x_test.shape[0], size=700, replace=False)

res.iloc[sample_indices,0:30].to_csv("/oak/stanford/groups/ebutcher/kevin/shap_test3_models.csv")

shap_values = explainer(x_test[sample_indices])

np.savetxt("/oak/stanford/groups/ebutcher/kevin/shap_test3.csv", np.squeeze(shap_values.values).T, delimiter=",")


shap.plots.beeswarm(shap_values)
















class SelfAttention(layers.Layer):
    def __init__(self, **kwargs):
        super(SelfAttention, self).__init__(**kwargs)
        self.supports_masking = True  # Allows handling variable-length sequences if applicable
    def build(self, input_shape):
        self.query_dense = layers.Dense(input_shape[-1], name="query_dense")
        self.key_dense = layers.Dense(input_shape[-1], name="key_dense")
        self.value_dense = layers.Dense(input_shape[-1], name="value_dense")
        self.combine_heads = layers.Dense(input_shape[-1], name="combine_heads")
        super(SelfAttention, self).build(input_shape)
    def call(self, inputs):
        # inputs shape: (batch_size, num_features, embedding_dim) or (batch_size, num_features) if no embeddings
        query = self.query_dense(inputs)
        key = self.key_dense(inputs)
        value = self.value_dense(inputs)
        # Scaled dot-product attention
        attention_scores = tf.matmul(query, key, transpose_b=True)
        dk = tf.cast(tf.shape(key)[-1], tf.float32)
        scaled_attention_scores = attention_scores / tf.math.sqrt(dk)
        attention_weights = tf.nn.softmax(scaled_attention_scores, axis=-1)
        # Apply attention weights to values
        output = tf.matmul(attention_weights, value)
        # Combine the outputs (this is for single-head attention, for multi-head, concatenation would happen here)
        output = self.combine_heads(output)
        return output
    def get_config(self):
        config = super().get_config()
        return config




numerical_inputs = keras.Input(shape=(len(numerical_features),), name="numerical_inputs")
categorical_inputs_1 = keras.Input(shape=(1,), name="categorical_inputs_1")
categorical_inputs_2 = keras.Input(shape=(1,), name="categorical_inputs_2")

# Embedding for categorical features
embedding_dim = 32  # Choose a suitable embedding dimension
embed_cat1 = layers.Embedding(input_dim=df['categorical_feature_1'].nunique(), output_dim=embedding_dim)(categorical_inputs_1)
embed_cat2 = layers.Embedding(input_dim=df['categorical_feature_2'].nunique(), output_dim=embedding_dim)(categorical_inputs_2)

# Flatten embeddings to combine with numerical features
embed_cat1_flat = layers.Flatten()(embed_cat1)
embed_cat2_flat = layers.Flatten()(embed_cat2)

# Concatenate all features
all_features = layers.concatenate([numerical_inputs, embed_cat1_flat, embed_cat2_flat])

# Reshape for self-attention (treating each feature as a "token")
# Assuming we want to attend over features, reshape to (batch_size, num_features, feature_dim)
# We can conceptualize each feature (numerical or embedded categorical) as a "token".
# For simplicity, we'll embed numerical features to the same dimension and concatenate,
# but a more sophisticated approach would involve feature tokenizers.
feature_dim = embedding_dim  # Assuming all feature representations will have this dimension

# Need to project numerical features to the same embedding_dim for concatenation
numerical_projected = layers.Dense(feature_dim)(numerical_inputs)

# Combine embeddings and projected numerical features (requires careful structuring)
# For this simple example, we'll just use the concatenated features directly with self-attention.
# In a true TabTransformer, embeddings would be fed into the Transformer layers and concatenated with numerical features later.

# Reshape for self-attention
# The `all_features` tensor is (batch_size, total_features_dimension).
# To apply self-attention over features, we need to treat each feature as an element in a sequence.
# Let's imagine each feature is represented by a vector.
# This requires a more complex input structure or a different application of self-attention.

# For a simplified approach, we'll assume `all_features` represents our "sequence" of features.
# This is a simplification and not exactly how TabTransformer processes inputs.
# A more accurate approach would involve transforming each categorical feature into an embedding,
# and then applying self-attention over these embeddings (potentially concatenated with numerical features after projection).

# Let's create a more TabTransformer-like structure:
# 1. Embed categorical features (already done)
# 2. Optionally, project numerical features to the same embedding dimension
# 3. Concatenate all feature embeddings/projections to create a sequence for self-attention.

# Let's rebuild the input processing for clarity in the context of TabTransformer:

# Numerical features
numerical_features_input = keras.Input(shape=(len(numerical_features),), name="numerical_features_input")

# Categorical feature embeddings
categorical_feature_1_input = keras.Input(shape=(1,), name="categorical_feature_1_input")
categorical_feature_2_input = keras.Input(shape=(1,), name="categorical_feature_2_input")

cat_embed_1 = layers.Embedding(input_dim=df['categorical_feature_1'].nunique(), output_dim=embedding_dim)(categorical_feature_1_input)
cat_embed_2 = layers.Embedding(input_dim=df['categorical_feature_2'].nunique(), output_dim=embedding_dim)(categorical_feature_2_input)

# Column embeddings (as in TabTransformer)
col_embed_1 = layers.Embedding(input_dim=1, output_dim=embedding_dim)(tf.constant([0])) # Dummy for first categorical feature
col_embed_2 = layers.Embedding(input_dim=1, output_dim=embedding_dim)(tf.constant([0])) # Dummy for second categorical feature

cat_embed_1_plus_col = layers.Add()([cat_embed_1, col_embed_1])
cat_embed_2_plus_col = layers.Add()([cat_embed_2, col_embed_2])

# Concatenate feature embeddings for the self-attention layer
# Stack these embeddings to form a "sequence" of features
# Shape will be (batch_size, num_categorical_features, embedding_dim)
stacked_categorical_embeddings = layers.concatenate([cat_embed_1_plus_col, cat_embed_2_plus_col], axis=1)

# Apply self-attention
attention_output = SelfAttention()(stacked_categorical_embeddings)

# Flatten the attention output and concatenate with numerical features
attention_output_flat = layers.Flatten()(attention_output)
final_features = layers.concatenate([attention_output_flat, numerical_features_input])

# MLP for final classification
mlp_output = layers.Dense(64, activation="relu")(final_features)
mlp_output = layers.Dropout(0.3)(mlp_output)
mlp_output = layers.Dense(32, activation="relu")(mlp_output)
output_layer = layers.Dense(1, activation="sigmoid")(mlp_output)

# Create the model
model = keras.Model(inputs=[numerical_features_input, categorical_feature_1_input, categorical_feature_2_input], outputs=output_layer)

# Compile and train
model.compile(optimizer="adam", loss="binary_crossentropy", metrics=["accuracy"])

# Prepare inputs for training (adjusting for the new input structure)
X_train_numerical = X_train[numerical_features].values
X_test_numerical = X_test[numerical_features].values

X_train_cat1 = X_train['categorical_feature_1'].values
X_test_cat1 = X_test['categorical_feature_1'].values

X_train_cat2 = X_train['categorical_feature_2'].values
X_test_cat2 = X_test['categorical_feature_2'].values

model.fit(
    [X_train_numerical, X_train_cat1, X_train_cat2], y_train,
    epochs=10,
    batch_size=32,
    validation_data=([X_test_numerical, X_test_cat1, X_test_cat2], y_test)
)


























import torch
import torch.nn as nn
from tab_transformer_pytorch import TabTransformer

cont_mean_std = torch.randn(10, 2)

model = TabTransformer(
    categories = (10, 5, 6, 5, 8),      # tuple containing the number of unique values within each category
    num_continuous = 10,                # number of continuous values
    dim = 32,                           # dimension, paper set at 32
    dim_out = 1,                        # binary prediction, but could be anything
    depth = 6,                          # depth, paper recommended 6
    heads = 8,                          # heads, paper recommends 8
    attn_dropout = 0.1,                 # post-attention dropout
    ff_dropout = 0.1,                   # feed forward dropout
    mlp_hidden_mults = (4, 2),          # relative multiples of each hidden dimension of the last mlp to logits
    mlp_act = nn.ReLU(),                # activation for final mlp, defaults to relu, but could be anything else (selu etc)
    continuous_mean_std = cont_mean_std # (optional) - normalize the continuous values before layer norm
)

x_categ = torch.randint(0, 5, (1, 5))     # category values, from 0 - max number of categories, in the order as passed into the constructor above
x_cont = torch.randn(1, 10)               # assume continuous values are already normalized individually

pred = model(x_categ, x_cont) # (1, 1)



