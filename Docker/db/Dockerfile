FROM mysql:5.7


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
    mysql -u root -p${MYSQL_ROOT_PASSWORD} -e "GRANT ALL PRIVILEGES ON *.* TO 'edivapublic'@'%'; ;" && \
    mysql -u root -p${MYSQL_ROOT_PASSWORD} -e "CREATE DATABASE eDiVa_public_omics " && \
    mysql -u root -p${MYSQL_ROOT_PASSWORD} -e "CREATE DATABASE eDiVa_annotation " 
VOLUME /var/lib/mysql2/
CMD ["true"]

