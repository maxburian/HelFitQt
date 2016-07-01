#ifndef DIALOG_MODELVARS_H
#define DIALOG_MODELVARS_H

#include <helfitqt.h>
#include "ui_dialog_modelvars.h"

class dialog_modelvars : public QDialog
{
	Q_OBJECT

public:
	dialog_modelvars(QWidget *parent = 0);
	~dialog_modelvars();
	bool dialogOK=false;
	int m_nr_atoms;
	double m_stackspacing;
	double m_diameter;

private slots:
	void on_okButton_clicked();
	void on_cancelButton_clicked();

private:
	Ui::dialog_modelvars_class *ui;
};

#endif // DIALOG_MODELVARS_H