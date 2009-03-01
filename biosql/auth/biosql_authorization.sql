-- MySQL database definitions for authorization.
-- This provides a simple system of users and groups that organize a set
-- of users. Group associations are then associated with a permission set
-- which can have one or more groups. 
--
-- This permission mechanism provides extra flexibility and prevents the
-- generation of groups just to handle a permission issue. For instance,
-- you might want to allow users in two labs to access a set of data without
-- creating another group to which they both belong.
--
-- These permissions can be checked either in code or via explicit database
-- associations.
-- 
-- This is based on the TurboGears python authorization framework:
-- http://turbogears.org/2.0/docs/main/Extensions/Authorization/index.html

CREATE TABLE auth_user(
  user_id INT(10) UNSIGNED NOT NULL auto_increment,
  user_name VARCHAR(64) NOT NULL,
  password VARCHAR(40) NOT NULL,
  email_address VARCHAR(255),
  display_name VARCHAR(255),
  PRIMARY KEY (user_id),
  UNIQUE (email_address),
  UNIQUE (user_name)
) TYPE=INNODB;

-- Provide key/value pairs that can be attached to a user. This is useful
-- for associating data like preferences or saved information with a particular
-- user of your system.
CREATE TABLE auth_user_qualifier_value (
    user_id int(10) unsigned default NULL,
    term_id       INT(10) UNSIGNED NOT NULL,
    value         TEXT,
    rank	  INT(5) NOT NULL DEFAULT 0,
    UNIQUE (user_id,term_id,rank),
    FOREIGN KEY (user_id) REFERENCES auth_user(user_id),
    FOREIGN KEY (term_id) REFERENCES term (term_id),
    KEY authuserqual_trm (term_id)
) TYPE=INNODB;

CREATE TABLE auth_group(
  group_id INT(10) UNSIGNED NOT NULL auto_increment,
  group_name VARCHAR(64) NOT NULL,
  description VARCHAR(255),
  PRIMARY KEY (group_id),
  UNIQUE (group_name)
) TYPE=INNODB;

CREATE TABLE auth_permission(
  permission_id INT(10) UNSIGNED NOT NULL auto_increment,
  permission_name VARCHAR(64) NOT NULL,
  description VARCHAR(255),
  PRIMARY KEY (permission_id),
  UNIQUE (permission_name)
) TYPE=INNODB;

CREATE TABLE auth_user_group(
  user_id INT(10) UNSIGNED NOT NULL,
  group_id INT(10) UNSIGNED NOT NULL,
  FOREIGN KEY (user_id) REFERENCES auth_user(user_id),
  FOREIGN KEY (group_id) REFERENCES auth_group(group_id)
) TYPE=INNODB;

CREATE TABLE auth_group_permission(
  group_id INT(10) UNSIGNED NOT NULL,
  permission_id INT(10) UNSIGNED NOT NULL,
  FOREIGN KEY (group_id) REFERENCES auth_group(group_id),
  FOREIGN KEY (permission_id) REFERENCES auth_permission(permission_id)
) TYPE=INNODB;
