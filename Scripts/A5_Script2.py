# BINF6210 Software Tools
# Assignment 5
# Script File 2 - Machine Learning Model Comparison for Biomarker Validation in Python
# Isha Baxi
# Optimizing Gene Expression Analysis for Alzheimer's Disease with Machine Learning for Biomarker Identification and Validation

# Import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier  # Random Forest classifier
from sklearn.svm import SVC  # Support Vector Classifier
from sklearn.model_selection import train_test_split  # For splitting dataset
from sklearn.preprocessing import StandardScaler, LabelEncoder  # For data preprocessing
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score, roc_curve  # Evaluation metrics
from imblearn.over_sampling import SMOTE  # To handle class imbalance through oversampling

# Load the gene expression data where the last column is the target 'Condition'
file_path = '../Data/ml_input_common_genes.csv'
data = pd.read_csv(file_path, index_col=0)

# Preprocessing data
X = data.iloc[:, :-1]  # Features (all columns except the target)
y = data['Condition']  # Target variable ('Condition')

# Encode the target variable (Condition) into numeric labels (e.g., 0 and 1)
label_encoder = LabelEncoder()
y_encoded = label_encoder.fit_transform(y)

# Split data into training (70%) and testing (30%) sets, ensuring stratification based on the target variable
X_train, X_test, y_train, y_test = train_test_split(
    X, y_encoded, test_size=0.3, random_state=42, stratify=y_encoded
)

# Handle class imbalance using SMOTE - SMOTE (Synthetic Minority Over-sampling Technique) is applied to the training data to balance the class distribution
smote = SMOTE(random_state=42)
X_train_resampled, y_train_resampled = smote.fit_resample(X_train, y_train)

# Standardize features to have a mean of 0 and a standard deviation of 1 for both training and test sets
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train_resampled)  # Fit on train, transform both train and test
X_test_scaled = scaler.transform(X_test)

# The following parameters were chosen based on suboptimal performance observed with the default model settings and insights drawn from relevant studies.

# Train Random Forest Classifier with specified parameters (balanced class weights)
# RF model with 5 estimators is set to reduce computational cost while giving enough tree diversity for a robust prediction. 
# The depth was set to 10 to limit growth of the decision tree to avoid overfitting, specifically because of the high-dim of gene expression data. 
# Class weights are balanced to manage class imbalance in sample, and random state is 42 to ensure same reproducibility across all runs. 
rf_model = RandomForestClassifier(
    n_estimators=5, max_depth=10, class_weight='balanced', random_state=42
)
rf_model.fit(X_train_scaled, y_train_resampled)

# Train SVM Model with linear kernel to classify the data
# Linear kernel is chosen for the SVM model to efficiently handle high-dimensional gene expression data while offering interpretability through feature coefficients. 

svm_model = SVC(
    kernel='linear',          # Linear kernel for SVM
    C = 1.0,                  # Regularization parameter, balances fitting training and avoiding overfitting
    class_weight='balanced',  # Handle class imbalance
    probability=True,         # Enable probability estimation for ROC-AUC calculation
    random_state=42           # Ensures same reproducibility across runs
)
svm_model.fit(X_train_scaled, y_train_resampled)

# Predictions and Metrics for Random Forest and SVM
# Make predictions and calculate evaluation metrics (accuracy, F1 score, AUC) for the models
y_pred_rf = rf_model.predict(X_test_scaled)  # Predicted labels
y_pred_proba_rf = rf_model.predict_proba(X_test_scaled)[:, 1]  # Predicted probabilities for ROC curve

rf_accuracy = accuracy_score(y_test, y_pred_rf)  # Accuracy score
rf_f1 = f1_score(y_test, y_pred_rf)  # F1 score
rf_auc = roc_auc_score(y_test, y_pred_proba_rf)  # AUC score

y_pred_svm = svm_model.predict(X_test_scaled)
y_pred_proba_svm = svm_model.predict_proba(X_test_scaled)[:, 1]

svm_accuracy = accuracy_score(y_test, y_pred_svm)
svm_f1 = f1_score(y_test, y_pred_svm)
svm_auc = roc_auc_score(y_test, y_pred_proba_svm)

# Feature Importance for Random Forest and SVM
# Get the feature importances from the Random Forest model (which can be used to interpret the model)
rf_feature_importance = pd.DataFrame({
    'Gene': X.columns,  # Genes are the features in the dataset
    'Importance': rf_model.feature_importances_
}).sort_values(by='Importance', ascending=False).head(10)  # Top 10 features with highest importance

# Get the feature importance for SVM model using the absolute value of the coefficients (for a linear kernel)
svm_feature_importance = pd.DataFrame({
    'Gene': X.columns,
    'Importance': np.abs(svm_model.coef_[0])  # SVM coefficients represent feature weights for linear kernel
}).sort_values(by='Importance', ascending=False).head(10)

# Display Results
# Print evaluation results (accuracy, F1 score, AUC) and top 10 most important features for both models
print("--- Random Forest ---")
print(f"Accuracy: {rf_accuracy:.4f}")
print(f"F1 Score: {rf_f1:.4f}")
print(f"AUC: {rf_auc:.4f}")
print("Top 10 most important features:")
print(rf_feature_importance)

print("\n--- SVM ---")
print(f"Accuracy: {svm_accuracy:.4f}")
print(f"F1 Score: {svm_f1:.4f}")
print(f"AUC: {svm_auc:.4f}")
print("Top 10 most important features:")
print(svm_feature_importance)

# Plot Side-by-Side Comparison of Top 10 Features
# In the figure, the Random Forest model, show genes such as ANXA6, USP17L2, and ATP5O with the highest feature importance, suggesting their strong predictive contribution. The SVM model emphasizes a different set of genes, with IRX2, NLRP4, and MGC34796 showing the greatest importance. The contrast between the two panels indicates that the two algorithms may prioritize distinct sets of genes based on their decision-making processes. Both need to be validated with studies to cross-reference which model may be the best predictor. 
plt.figure(figsize=(14, 8))
sns.set_style("whitegrid")

# Barplot for Random Forest 
plt.subplot(1, 2, 1)
sns.barplot(
    data=rf_feature_importance,
    y='Gene', x='Importance',
    palette='Greens_r', edgecolor="black"
)
plt.title("Top 10 Genes - Random Forest", fontsize=14)
plt.xlabel("Feature Importance", fontsize=12)
plt.ylabel("Gene", fontsize=12)
plt.gca().invert_yaxis()  # Highest importance on top

# Barplot for SVM 
plt.subplot(1, 2, 2)
sns.barplot(
    data=svm_feature_importance,
    y='Gene', x='Importance',
    palette='Blues_r', edgecolor="black"
)
plt.title("Top 10 Genes - SVM", fontsize=14)
plt.xlabel("Feature Importance", fontsize=12)
plt.ylabel("")  # Remove redundant ylabel
plt.gca().invert_yaxis()  # Highest importance on top

plt.tight_layout()  # Ensure proper spacing between plots
plt.show()

# Plot ROC Curves for both models, generate and plot ROC curves to assess the performance of both models in terms of true positive rate and false positive rate
# In the figure, the SVM model (blue curve) demonstrates superior performance with an AUC of 0.90, indicating strong discriminatory ability. In contrast, the Random Forest model (green curve) has a lower AUC of 0.70, suggesting relatively weaker predictive performance. This comparison highlights the effectiveness of SVM in identifying Alzheimer’s-related patterns in the gene expression data.

# However, when we consider the overall performance metrics, additional nuances emerge.
# For the RF the accuracy is performing slightly better than SVM, but it is only slightly better than random chance. When looking at F1 score, again RF performs slightly better indicating decent balance between precision and recall, suggesting it is reasonably good at correctly identifying AD-related cases while minimizing false positives and false negatives.
# These metrics reveal a trade-off between the models. While SVM has higher AUC and focuses on distinguishing between classes effectively, its low accuracy and F1 score suggest practical limitations in application. Alternatively, RF demonstrates more balanced performance metrics, with a reasonable F1 score and accuracy, but lacks the high discriminatory power of SVM as reflected by the AUC.

# Random Forest and SVM ROC Curve
fpr_rf, tpr_rf, _ = roc_curve(y_test, y_pred_proba_rf)
fpr_svm, tpr_svm, _ = roc_curve(y_test, y_pred_proba_svm)

# Plot both ROC curves on the same graph
plt.figure(figsize=(10, 7))
plt.plot(fpr_rf, tpr_rf, color='green', lw=2, label=f"Random Forest (AUC = {rf_auc:.2f})")
plt.plot(fpr_svm, tpr_svm, color='blue', lw=2, label=f"SVM (AUC = {svm_auc:.2f})")
plt.plot([0, 1], [0, 1], color='orange', lw=2, linestyle='--')  # Diagonal line (random classifier)
plt.xlabel("False Positive Rate", fontsize=14)
plt.ylabel("True Positive Rate", fontsize=14)
plt.title("ROC Curves for Alzheimer's Prediction Models", fontsize=16)
plt.legend(loc="lower right")  # Display legend in the lower-right corner
plt.grid(alpha=0.3)  # Add grid for better readability
plt.tight_layout()  # Ensure tight layout
plt.show()