/********************************************************************************
** Form generated from reading UI file 'dialog_modelvars.ui'
**
** Created by: Qt User Interface Compiler version 5.6.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_DIALOG_MODELVARS_H
#define UI_DIALOG_MODELVARS_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_dialog_modelvars_class
{
public:
    QWidget *layoutWidget;
    QVBoxLayout *vboxLayout;
    QSpacerItem *spacer;
    QPushButton *okButton;
    QPushButton *cancelButton;
    QSpacerItem *spacerItem;
    QSpinBox *spbNrAtoms;
    QLineEdit *lineStackSpacing;
    QLineEdit *lineDiameter;
    QLabel *label;
    QLabel *label_2;
    QLabel *label_3;
    QLabel *label_24;
    QLabel *label_25;
    QLabel *label_15;

    void setupUi(QDialog *dialog_modelvars_class)
    {
        if (dialog_modelvars_class->objectName().isEmpty())
            dialog_modelvars_class->setObjectName(QStringLiteral("dialog_modelvars_class"));
        dialog_modelvars_class->resize(400, 159);
        dialog_modelvars_class->setModal(true);
        layoutWidget = new QWidget(dialog_modelvars_class);
        layoutWidget->setObjectName(QStringLiteral("layoutWidget"));
        layoutWidget->setGeometry(QRect(290, 25, 95, 131));
        vboxLayout = new QVBoxLayout(layoutWidget);
        vboxLayout->setSpacing(6);
        vboxLayout->setObjectName(QStringLiteral("vboxLayout"));
        vboxLayout->setContentsMargins(0, 0, 0, 0);
        spacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        vboxLayout->addItem(spacer);

        okButton = new QPushButton(layoutWidget);
        okButton->setObjectName(QStringLiteral("okButton"));

        vboxLayout->addWidget(okButton);

        cancelButton = new QPushButton(layoutWidget);
        cancelButton->setObjectName(QStringLiteral("cancelButton"));

        vboxLayout->addWidget(cancelButton);

        spacerItem = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        vboxLayout->addItem(spacerItem);

        spbNrAtoms = new QSpinBox(dialog_modelvars_class);
        spbNrAtoms->setObjectName(QStringLiteral("spbNrAtoms"));
        spbNrAtoms->setGeometry(QRect(140, 50, 61, 22));
        spbNrAtoms->setMinimum(1);
        spbNrAtoms->setMaximum(5000);
        spbNrAtoms->setValue(500);
        lineStackSpacing = new QLineEdit(dialog_modelvars_class);
        lineStackSpacing->setObjectName(QStringLiteral("lineStackSpacing"));
        lineStackSpacing->setGeometry(QRect(140, 80, 61, 22));
        lineDiameter = new QLineEdit(dialog_modelvars_class);
        lineDiameter->setObjectName(QStringLiteral("lineDiameter"));
        lineDiameter->setGeometry(QRect(140, 110, 61, 22));
        label = new QLabel(dialog_modelvars_class);
        label->setObjectName(QStringLiteral("label"));
        label->setGeometry(QRect(10, 50, 121, 20));
        label->setLayoutDirection(Qt::LeftToRight);
        label->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_2 = new QLabel(dialog_modelvars_class);
        label_2->setObjectName(QStringLiteral("label_2"));
        label_2->setGeometry(QRect(10, 80, 121, 20));
        label_2->setLayoutDirection(Qt::LeftToRight);
        label_2->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_3 = new QLabel(dialog_modelvars_class);
        label_3->setObjectName(QStringLiteral("label_3"));
        label_3->setGeometry(QRect(10, 110, 121, 20));
        label_3->setLayoutDirection(Qt::LeftToRight);
        label_3->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_24 = new QLabel(dialog_modelvars_class);
        label_24->setObjectName(QStringLiteral("label_24"));
        label_24->setGeometry(QRect(210, 80, 51, 16));
        label_25 = new QLabel(dialog_modelvars_class);
        label_25->setObjectName(QStringLiteral("label_25"));
        label_25->setGeometry(QRect(210, 110, 51, 16));
        label_15 = new QLabel(dialog_modelvars_class);
        label_15->setObjectName(QStringLiteral("label_15"));
        label_15->setGeometry(QRect(30, 10, 311, 20));
        QFont font;
        font.setPointSize(10);
        font.setBold(true);
        font.setWeight(75);
        label_15->setFont(font);
#ifndef QT_NO_SHORTCUT
        label_24->setBuddy(lineStackSpacing);
        label_25->setBuddy(lineStackSpacing);
#endif // QT_NO_SHORTCUT

        retranslateUi(dialog_modelvars_class);

        QMetaObject::connectSlotsByName(dialog_modelvars_class);
    } // setupUi

    void retranslateUi(QDialog *dialog_modelvars_class)
    {
        dialog_modelvars_class->setWindowTitle(QApplication::translate("dialog_modelvars_class", "initial seed parameters...", 0));
        okButton->setText(QApplication::translate("dialog_modelvars_class", "OK", 0));
        cancelButton->setText(QApplication::translate("dialog_modelvars_class", "Cancel", 0));
        lineStackSpacing->setText(QApplication::translate("dialog_modelvars_class", "25.0", 0));
        lineDiameter->setText(QApplication::translate("dialog_modelvars_class", "20", 0));
        label->setText(QApplication::translate("dialog_modelvars_class", "Number of atoms:", 0));
        label_2->setText(QApplication::translate("dialog_modelvars_class", "height / stack dist.:", 0));
        label_3->setText(QApplication::translate("dialog_modelvars_class", "diameter:", 0));
        label_24->setText(QApplication::translate("dialog_modelvars_class", "(nm)", 0));
        label_25->setText(QApplication::translate("dialog_modelvars_class", "(nm)", 0));
        label_15->setText(QApplication::translate("dialog_modelvars_class", "Specify parameters for initial model..", 0));
    } // retranslateUi

};

namespace Ui {
    class dialog_modelvars_class: public Ui_dialog_modelvars_class {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_DIALOG_MODELVARS_H
