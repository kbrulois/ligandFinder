

import tensorflow as tf
from tensorflow.keras.layers import (
    Input, Dense, Dropout, LayerNormalization,
    MultiHeadAttention, Add, Flatten
)
from tensorflow.keras.models import Model

from tensorflow.keras.regularizers import l2

def transformer_block(inputs, num_heads=2, ff_dim=128, dropout_rate=0.5):
    attention_output = MultiHeadAttention(num_heads=num_heads, key_dim=inputs.shape[-1])(inputs, inputs)
    attention_output = Dropout(dropout_rate)(attention_output)
    out1 = LayerNormalization(epsilon=1e-6)(Add()([inputs, attention_output]))
    ffn_output = Dense(ff_dim, activation='relu')(out1)
    ffn_output = Dense(inputs.shape[-1])(ffn_output)
    ffn_output = Dropout(dropout_rate)(ffn_output)
    out2 = LayerNormalization(epsilon=1e-6)(Add()([out1, ffn_output]))
    return out2
  
  
def build_transformer_model(input_dim, num_heads=2, ff_dim=128):
    inputs = Input(shape=(input_dim,))
    x = Dense(ff_dim, activation='relu', kernel_regularizer=l2(1e-4))(inputs)
    x = tf.expand_dims(x, axis=1)  # (batch_size, sequence_len=1, embedding_dim)
    x = transformer_block(x, num_heads=num_heads, ff_dim=ff_dim)
    x = Flatten()(x)
    x = Dense(32, activation='relu', kernel_regularizer=l2(1e-4))(x)
    x = Dropout(0.2)(x)
    outputs = Dense(1, activation='sigmoid')(x)
    return Model(inputs=inputs, outputs=outputs)
  
  def build_attention_model_no_projection(input_dim, num_heads=2, ff_dim=64):
    inputs = Input(shape=(input_dim,))  # flat input, e.g., (35,)
    x = tf.expand_dims(inputs, axis=-1)  # shape becomes (batch, 35, 1)
    x = transformer_block(x, num_heads=num_heads, ff_dim=ff_dim)
    x = Flatten()(x)
    x = Dense(32, activation='relu')(x)
    x = Dropout(0.1)(x)
    outputs = Dense(1, activation='sigmoid')(x)
    return Model(inputs, outputs)

# Example usage:
model = build_transformer_model(input_dim=x_train.shape[1])

model.compile(optimizer=keras.optimizers.Adam(learning_rate=1e-1),
              loss='binary_crossentropy',
              metrics=metrics)


model.fit(x_train, 
          y_train, 
          batch_size = 64,
          epochs=100,
          callbacks=callbacks,
          validation_data=(x_val, y_val),
          class_weight=class_weight)


model.evaluate(x_test, y_test, verbose=2)
