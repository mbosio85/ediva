


FROM mysql:5.7

#EXPOSE 3306
MAINTAINER mattia bosio

# Environment variables
ENV MYSQL_ROOT_PASSWORD=rootpwd
ENV MYSQL_HOST=localhost
# Allows you to change the value of "max_allowed_packet"
# ["mysqlconf/gatewaymy.cnf", "/etc/mysql/conf.d/conf_mysql.cnf"]

# Create Database
RUN	mkdir /usr/sql
RUN	chmod 644 /usr/sql

RUN sed -i 's|/var/lib/mysql|/var/lib/mysql2|g' /etc/mysql/mysql.conf.d/mysqld.cnf

RUN echo $MYSQL_ROOT_PASSWORD 
#RUN sed -i 's|/var/lib/mysql|/var/lib/mysql2|g' /etc/etc/mysql/my.cnf

#add user edivapublib with password:
RUN /etc/init.d/mysql start && \
    mysql -u root -p${MYSQL_ROOT_PASSWORD} -e "CREATE USER 'edivapublic'@'%' IDENTIFIED BY 'x86d2k1B';" && \
    mysql -u root -p${MYSQL_ROOT_PASSWORD} -e "GRANT ALL PRIVILEGES ON *.* TO 'edivapublic'@'%' ;" && \
    mysql -u root -p${MYSQL_ROOT_PASSWORD} -e "CREATE DATABASE eDiVa_public_omics; " && \
    mysql -u root -p${MYSQL_ROOT_PASSWORD} -e "CREATE DATABASE eDiVa_annotation; " 


# This does not work because it exceeds Docker image size limits.
# So we came up with a data volume to overcome this limitation.
# If you can alter your default data you can directly build the image from the Dockerfile.
#eDiVA_public_omics
#COPY ["eDiVa_public_omics.sql.gz", "/usr/sql/eDiVa_public_omics.sql.gz"]
#RUN /etc/init.d/mysql start && \
#    	zcat /usr/sql/eDiVa_public_omics.sql.gz | mysql -u edivapublic  -px86d2k1B -D eDiVa_public_omics && \
#        rm -rd /usr/sql 
#eDiVA_annotation
# THIS IS GOING TO TAKE A LONG TIME #
#COPY ["eDiVa_annotation.sql.gz","/usr/sql/eDiVa_annotation.sql.gz"]
#RUN     echo "THIS IS GOING TO TAKE A LONG TIME LOADING ANNOTATION DATABASE"
#        /etc/init.d/mysql start && \
#        zcat /usr/sql/eDiVa_annotation.sql.gz mysql -u edivapublic -px86d2k1B -D eDiVa_annotation && \
#        rm -rd /usr/sql 




